"""Omics integration and energy analysis (Phase 5).

Integrate gene-expression (transcriptomics/proteomics) data with the model to
build *context-specific* predictions, plus an ATP-maintenance sensitivity scan.

* **eFlux** [Colijn et al. 2009] — scales each reaction's flux bounds by the
  expression of the enzymes that catalyze it (via the gene–reaction rule), so the
  flux capacity reflects how strongly the gene set is expressed.
* **GIMME** [Becker & Palsson 2008] — finds a flux distribution that satisfies a
  required objective while minimizing flux through reactions whose genes are
  expressed below a threshold (an inconsistency score quantifies the conflict).
* **ATP-maintenance sensitivity** — scans the non-growth ATP maintenance demand
  and reports the effect on growth, revealing the energetic burden a strain can
  tolerate.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra
import numpy as np
import pandas as pd
from optlang.symbolics import Zero

from ..editing import guess_biomass_reaction


class OmicsError(Exception):
    """Raised when omics integration cannot be performed."""


@dataclass
class OmicsResult:
    method: str
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    note: str = ""
    objective_value: float = float("nan")


# --- expression data --------------------------------------------------------
def load_expression(path: str) -> Dict[str, float]:
    """Load a two-column gene-expression table (gene id, value) from CSV/TSV."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception as exc:  # noqa: BLE001
        raise OmicsError(f"Could not read expression file:\n{exc}") from exc
    if df.shape[1] < 2:
        raise OmicsError("Expression file needs at least two columns: gene id and value.")
    genes = df.iloc[:, 0].astype(str)
    try:
        values = df.iloc[:, 1].astype(float)
    except Exception as exc:  # noqa: BLE001
        raise OmicsError(f"The second column must be numeric expression values:\n{exc}") from exc
    return dict(zip(genes, values))


def reaction_expression(model: cobra.Model, expression: Dict[str, float],
                        default: Optional[float] = None) -> Dict[str, float]:
    """Map gene expression to a value per reaction via the gene-reaction rule.

    AND (enzyme complex) -> minimum of the components; OR (isozymes) -> sum.
    Reactions without a usable rule get ``default`` (the median expression if None).
    """
    if default is None:
        default = float(np.median(list(expression.values()))) if expression else 0.0
    result = {}
    for rxn in model.reactions:
        rule = rxn.gene_reaction_rule.strip()
        if not rule:
            result[rxn.id] = default
            continue
        result[rxn.id] = _eval_gpr(rule, expression, default)
    return result


def _eval_gpr(rule: str, expression: Dict[str, float], default: float) -> float:
    """Evaluate a GPR string to a number (AND=min, OR=sum)."""
    # Normalize boolean operators and tokenize.
    expr = rule.replace(" AND ", " and ").replace(" OR ", " or ")

    def repl(match):
        gid = match.group(0)
        if gid in ("and", "or", "(", ")"):
            return gid
        return repr(float(expression.get(gid, default)))

    tokens = re.sub(r"[A-Za-z0-9_.\-:]+", repl, expr)
    # Convert to nested min()/sum() — evaluate safely via a tiny recursive parser.
    try:
        return _safe_eval_boolean(tokens)
    except Exception:  # noqa: BLE001
        return default


def _safe_eval_boolean(expr: str) -> float:
    """Evaluate an expression of numbers combined with 'and'/'or' and parentheses.

    'and' -> min, 'or' -> sum. Implemented with Python's parser on a restricted AST.
    """
    import ast

    node = ast.parse(expr, mode="eval").body

    def ev(n):
        if isinstance(n, ast.BoolOp):
            vals = [ev(v) for v in n.values]
            return min(vals) if isinstance(n.op, ast.And) else sum(vals)
        if isinstance(n, ast.Constant):
            return float(n.value)
        if isinstance(n, ast.UnaryOp):  # e.g. negative numbers
            return -ev(n.operand) if isinstance(n.op, ast.USub) else ev(n.operand)
        raise ValueError("unsupported node")

    return float(ev(node))


# --- eFlux ------------------------------------------------------------------
def run_eflux(model: cobra.Model, expression: Dict[str, float]) -> OmicsResult:
    """Constrain reaction bounds by expression and optimize the current objective."""
    rxn_expr = reaction_expression(model, expression)
    max_expr = max(rxn_expr.values()) if rxn_expr else 1.0
    if max_expr <= 0:
        raise OmicsError("All reaction expression values are zero.")

    work = model.copy()
    for rxn in work.reactions:
        if rxn.boundary:
            continue  # leave exchange/medium bounds untouched
        scale = rxn_expr.get(rxn.id, 0.0) / max_expr
        ub = abs(rxn.upper_bound) * scale
        lb = -abs(rxn.lower_bound) * scale if rxn.lower_bound < 0 else 0.0
        rxn.upper_bound = max(ub, 0.0)
        rxn.lower_bound = min(lb, 0.0) if rxn.lower_bound < 0 else 0.0

    sol = work.optimize()
    rows = [{"reaction": rid, "expression": rxn_expr.get(rid, 0.0),
             "flux": float(sol.fluxes.get(rid, 0.0)) if sol.status == "optimal" else float("nan")}
            for rid in (r.id for r in model.reactions)]
    table = pd.DataFrame(rows)
    obj = float(sol.objective_value) if sol.status == "optimal" else float("nan")
    return OmicsResult(method="eFlux", table=table,
                       objective_value=obj,
                       note=f"Status: {sol.status}; bounds scaled by gene expression.")


# --- GIMME ------------------------------------------------------------------
def run_gimme(model: cobra.Model, expression: Dict[str, float], *,
              threshold: float, objective_fraction: float = 0.9) -> OmicsResult:
    """GIMME: satisfy the objective while minimizing low-expression flux usage."""
    rxn_expr = reaction_expression(model, expression)
    biomass = guess_biomass_reaction(model)
    if not biomass:
        raise OmicsError("Could not identify the objective/biomass reaction.")

    work = model.copy()
    # Required metabolic functionality: objective >= fraction * optimum.
    opt = work.slim_optimize()
    if opt is None or not np.isfinite(opt):
        raise OmicsError("The model is infeasible under the current medium.")
    obj_rxn = work.reactions.get_by_id(biomass)
    work.problem  # ensure solver built
    rmf = work.problem.Constraint(obj_rxn.flux_expression,
                                  lb=objective_fraction * opt, name="gimme_rmf")
    work.add_cons_vars([rmf])

    # Minimize sum of penalty * |flux| for below-threshold reactions.
    coefficients = {}
    for rxn in work.reactions:
        penalty = threshold - rxn_expr.get(rxn.id, threshold)
        if penalty > 0:
            coefficients[rxn.forward_variable] = penalty
            coefficients[rxn.reverse_variable] = penalty
    if not coefficients:
        raise OmicsError("No reactions fall below the expression threshold; "
                         "lower the threshold to get a GIMME result.")
    objective = work.problem.Objective(Zero, direction="min")
    work.objective = objective
    work.objective.set_linear_coefficients(coefficients)

    sol = work.optimize()
    if sol.status != "optimal":
        raise OmicsError(f"GIMME optimization was not optimal (status: {sol.status}).")
    inconsistency = float(sol.objective_value)

    rows = [{"reaction": rid,
             "expression": rxn_expr.get(rid, 0.0),
             "below_threshold": rxn_expr.get(rid, threshold) < threshold,
             "flux": float(sol.fluxes.get(rid, 0.0))}
            for rid in (r.id for r in model.reactions)]
    table = pd.DataFrame(rows)
    return OmicsResult(method="GIMME", table=table, objective_value=inconsistency,
                       note=(f"Inconsistency score = {inconsistency:.4g} (lower is better); "
                             f"objective held ≥ {objective_fraction:.0%} of optimum."))


# --- ATP maintenance sensitivity --------------------------------------------
def run_atpm_sensitivity(model: cobra.Model, *, atpm_id: str = "ATPM",
                         points: int = 20, max_atpm: Optional[float] = None) -> OmicsResult:
    """Scan the ATP maintenance lower bound and record growth at each level."""
    if not model.reactions.has_id(atpm_id):
        raise OmicsError(f"No ATP-maintenance reaction '{atpm_id}' in the model. "
                         "Specify the correct reaction id.")
    work = model.copy()
    atpm = work.reactions.get_by_id(atpm_id)
    if max_atpm is None:
        # Find the largest ATPM the model can sustain.
        with work:
            work.objective = atpm
            top = work.slim_optimize()
        max_atpm = float(top) if top and np.isfinite(top) else 100.0
    values = np.linspace(0.0, max_atpm, points)
    rows = []
    for v in values:
        with work:
            atpm.lower_bound = v
            g = work.slim_optimize()
            rows.append({"atp_maintenance": float(v),
                         "growth": float(g) if g is not None and np.isfinite(g) else 0.0})
    return OmicsResult(method="ATPM sensitivity", table=pd.DataFrame(rows),
                       note=f"Growth vs non-growth ATP maintenance ({atpm_id}).")
