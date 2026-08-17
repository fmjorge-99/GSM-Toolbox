"""Model quality control — lightweight, built-in curation checks.

Inspired by the kinds of tests MEMOTE runs, but implemented directly on top of
COBRApy so the app stays self-contained (no heavy extra dependency). Reports:

* mass/charge-imbalanced reactions (excluding boundary/exchange/biomass),
* reactions whose balance **cannot be checked** because a participant has no formula,
* reactions with no gene-reaction rule (GPR),
* dead-end / orphan metabolites (only produced or only consumed),
* blocked reactions (cannot carry any flux), optional as it is slower.

Returns clean DataFrames plus a summary the GUI renders as a scorecard.

The first two checks are deliberately separate. cobra's ``check_mass_balance`` counts a
formula-less metabolite as contributing no atoms, so a reaction whose product has no
formula is reported as losing that product's entire mass — and the score is docked for a
fault that does not exist. Roughly 40% of the merged database's metabolites have no
formula, so the distinction is the difference between a usable report and a wall of
false positives.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import cobra
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions


@dataclass
class QCReport:
    summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    imbalanced: pd.DataFrame = field(default_factory=pd.DataFrame)
    #: reactions that could not be checked at all — a participant has no formula.
    #: Kept apart from `imbalanced`: this is a data gap, not a broken reaction.
    unverifiable: pd.DataFrame = field(default_factory=pd.DataFrame)
    missing_gpr: pd.DataFrame = field(default_factory=pd.DataFrame)
    dead_ends: pd.DataFrame = field(default_factory=pd.DataFrame)
    blocked: pd.DataFrame = field(default_factory=pd.DataFrame)
    score: float = 0.0


def _balance_checks(model: cobra.Model):
    """``(imbalanced, unverifiable)`` DataFrames for the model's metabolic reactions."""
    from ..balancing import assess_reaction

    bad, unknown = [], []
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        if "biomass" in rxn.id.lower():
            continue  # biomass is intentionally not mass-balanced
        try:
            verdict = assess_reaction(rxn)
        except Exception:  # noqa: BLE001 - malformed formula (R groups, polymers)
            unknown.append({"reaction": rxn.id,
                            "reason": "its formulas could not be parsed"})
            continue
        if verdict.verifiably_imbalanced:
            bad.append({
                "reaction": rxn.id,
                "imbalance": ", ".join(f"{k}: {v:g}" for k, v in
                                       sorted(verdict.residual.items())),
            })
        elif not verdict.checkable:
            unknown.append({"reaction": rxn.id,
                            "reason": verdict.missing_formula_sentence(max_named=4)})
    return pd.DataFrame(bad), pd.DataFrame(unknown)


def _missing_gpr(model: cobra.Model) -> pd.DataFrame:
    rows = [{"reaction": r.id, "name": r.name or ""}
            for r in model.reactions
            if not r.boundary and not r.gene_reaction_rule.strip()]
    return pd.DataFrame(rows)


def _dead_end_metabolites(model: cobra.Model) -> pd.DataFrame:
    """Metabolites that can only be produced or only be consumed (network gaps)."""
    rows = []
    for met in model.metabolites:
        produced = consumed = False
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            lo, hi = rxn.lower_bound, rxn.upper_bound
            # Account for reversibility when deciding production/consumption.
            if coeff > 0:
                produced |= hi > 0
                consumed |= lo < 0
            elif coeff < 0:
                consumed |= hi > 0
                produced |= lo < 0
        if not (produced and consumed):
            kind = "only produced" if produced else ("only consumed" if consumed else "isolated")
            rows.append({"metabolite": met.id, "name": met.name or "", "issue": kind})
    return pd.DataFrame(rows)


def quality_report(model: cobra.Model, *, include_blocked: bool = True) -> QCReport:
    """Run the QC checks and assemble a :class:`QCReport`."""
    imbalanced, unverifiable = _balance_checks(model)
    missing_gpr = _missing_gpr(model)
    dead_ends = _dead_end_metabolites(model)

    blocked = pd.DataFrame()
    if include_blocked:
        try:
            blocked_ids = find_blocked_reactions(model)
            blocked = pd.DataFrame({"reaction": list(blocked_ids)})
        except Exception:  # noqa: BLE001
            blocked = pd.DataFrame()

    n_rxn = max(len(model.reactions), 1)
    n_met = max(len(model.metabolites), 1)
    # Simple 0–100 score: fraction of clean reactions/metabolites across checks.
    # `unverifiable` is NOT penalised — docking the score for a formula the database
    # never shipped punishes the user for someone else's gap, and used to make a fine
    # model built on the merged database score near zero. It is reported instead.
    penalties = (
        0.4 * len(imbalanced) / n_rxn
        + 0.2 * len(missing_gpr) / n_rxn
        + 0.2 * len(dead_ends) / n_met
        + 0.2 * (len(blocked) / n_rxn if not blocked.empty else 0.0)
    )
    score = round(max(0.0, 1.0 - penalties) * 100, 1)

    summary = pd.DataFrame([
        {"check": "Overall quality score", "result": f"{score} / 100"},
        {"check": "Reactions", "result": len(model.reactions)},
        {"check": "Metabolites", "result": len(model.metabolites)},
        {"check": "Genes", "result": len(model.genes)},
        {"check": "Mass/charge-imbalanced reactions", "result": len(imbalanced)},
        {"check": "Balance not checkable (a metabolite has no formula)",
         "result": len(unverifiable)},
        {"check": "Reactions without GPR", "result": len(missing_gpr)},
        {"check": "Dead-end / orphan metabolites", "result": len(dead_ends)},
        {"check": "Blocked reactions", "result": len(blocked) if not blocked.empty else "not checked"},
    ])
    return QCReport(
        summary=summary, imbalanced=imbalanced, unverifiable=unverifiable,
        missing_gpr=missing_gpr, dead_ends=dead_ends, blocked=blocked, score=score,
    )
