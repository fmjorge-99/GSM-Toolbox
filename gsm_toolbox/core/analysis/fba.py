"""Flux Balance Analysis and its close relatives.

Phase 1 exposes FBA, parsimonious FBA (pFBA) and Flux Variability Analysis (FVA).
Results are returned as plain dataclasses / pandas objects so the GUI never has
to touch the solver directly.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List, Optional

import cobra
import pandas as pd
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis, pfba

# FVA solves two LPs per reaction. On Windows, spawning a worker pool (re-import
# + per-worker model pickling) is far costlier than the solves themselves, making
# cobra's default multiprocessing several times slower than one fast process.
# Force single-process on Windows; allow cobra's parallelism elsewhere.
_FVA_PROCESSES = 1 if os.name == "nt" else None


@dataclass
class FBAResult:
    """Result of an FBA / pFBA optimization."""

    status: str
    objective_value: float
    method: str = "FBA"
    fluxes: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    shadow_prices: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    reduced_costs: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    # Human-readable disclosure of any engineering compromise made to get a result
    # (e.g. pFBA had to relax the optimum). Empty when nothing notable happened.
    note: str = ""

    @property
    def is_optimal(self) -> bool:
        return self.status == "optimal"

    def active_fluxes(self, tol: float = 1e-9) -> pd.Series:
        """Fluxes whose magnitude exceeds ``tol`` (the reactions actually carrying flux)."""
        return self.fluxes[self.fluxes.abs() > tol].sort_values(key=abs, ascending=False)

    def flux_table(self) -> pd.DataFrame:
        """Return fluxes (and reduced costs when available) as a tidy DataFrame."""
        df = self.fluxes.rename("flux").to_frame()
        df.index.name = "reaction"
        if not self.reduced_costs.empty:
            df["reduced_cost"] = self.reduced_costs
        return df.reset_index()


def run_fba(model: cobra.Model) -> FBAResult:
    """Standard FBA. Maximizes (or minimizes) the current model objective."""
    solution = model.optimize()
    return _solution_to_result(solution, method="FBA", model=model)


def run_pfba(model: cobra.Model) -> FBAResult:
    """Parsimonious FBA: optimal objective with minimal total flux.

    Note: ``solution.objective_value`` from cobra's pFBA is the *minimized total
    flux*, not the growth rate. We report the model's true objective (computed
    from the returned fluxes) so the GUI shows growth, consistent with FBA.
    """
    # cobra's pfba() raises (rather than returning a status) when the model is
    # infeasible. On numerically poorly-scaled genome-scale models, requiring the
    # objective to *exactly* equal the optimum is often infeasible under GLPK's
    # tolerances even though FBA succeeds; a tiny relaxation fixes that. Fall back
    # gracefully (non-optimal result) if it is still infeasible.
    for fraction in (1.0, 0.9999):
        try:
            solution = pfba(model, fraction_of_optimum=fraction)
            result = _solution_to_result(solution, method="pFBA", model=model)
            if fraction < 1.0 and result.is_optimal:
                # Disclose the relaxation instead of silently accepting a sub-optimum.
                pct = (1.0 - fraction) * 100.0
                result.note = (
                    f"pFBA optimum was relaxed to {fraction:.4g}× the FBA optimum "
                    f"(within {pct:.2g}%): requiring the exact optimum was infeasible "
                    "under the solver's numerical tolerances on this model.")
            return result
        except OptimizationError:
            continue
    status = getattr(getattr(model, "solver", None), "status", "infeasible") or "infeasible"
    return FBAResult(status=str(status), objective_value=float("nan"), method="pFBA",
                     note="pFBA could not be solved even after relaxing the optimum.")


def run_fva(
    model: cobra.Model,
    *,
    fraction_of_optimum: float = 1.0,
    reaction_list: Optional[List[str]] = None,
    loopless: bool = False,
) -> pd.DataFrame:
    """Flux Variability Analysis.

    Returns a DataFrame indexed by reaction id with ``minimum``/``maximum`` columns.
    ``fraction_of_optimum`` constrains the objective to at least that fraction of
    its optimal value while exploring each reaction's flux range.
    """
    reactions = None
    if reaction_list:
        reactions = [model.reactions.get_by_id(r) for r in reaction_list]
    # cobra deprecated the boolean `loopless`; map to its accepted values.
    loopless_arg = "cycleFreeFlux" if loopless else None
    fva = flux_variability_analysis(
        model,
        reaction_list=reactions,
        fraction_of_optimum=fraction_of_optimum,
        loopless=loopless_arg,
        processes=_FVA_PROCESSES,
    )
    fva.index.name = "reaction"
    if _FVA_PROCESSES == 1:
        # Disclose the Windows single-process trade-off (see _FVA_PROCESSES note).
        fva.attrs["note"] = (
            "FVA ran single-process (Windows): here spawning worker processes "
            "re-imports COBRA and pickles the model per worker, which is slower "
            "than one fast process. This does not affect the results.")
    return fva


def _objective_from_fluxes(model: cobra.Model, fluxes: pd.Series) -> float:
    """Compute the model's objective value from a flux vector.

    Uses each reaction's ``objective_coefficient`` so the reported value is the
    real objective (e.g. growth), regardless of what the solver's internal
    objective was (pFBA minimizes total flux instead).
    """
    total = 0.0
    for rxn in model.reactions:
        coeff = rxn.objective_coefficient
        if coeff:
            total += coeff * float(fluxes.get(rxn.id, 0.0))
    return total


def _solution_to_result(solution: cobra.Solution, method: str, model: cobra.Model) -> FBAResult:
    if solution.status != "optimal":
        return FBAResult(status=solution.status, objective_value=float("nan"), method=method)
    # Shadow prices / reduced costs are not always populated (e.g. pFBA); guard them.
    try:
        shadow = solution.shadow_prices
    except Exception:  # noqa: BLE001
        shadow = pd.Series(dtype=float)
    try:
        reduced = solution.reduced_costs
    except Exception:  # noqa: BLE001
        reduced = pd.Series(dtype=float)
    objective_value = _objective_from_fluxes(model, solution.fluxes)
    return FBAResult(
        status=solution.status,
        objective_value=objective_value,
        method=method,
        fluxes=solution.fluxes,
        shadow_prices=shadow if shadow is not None else pd.Series(dtype=float),
        reduced_costs=reduced if reduced is not None else pd.Series(dtype=float),
    )
