"""Gap-filling: suggest reactions to add so the model can do something it can't.

Wraps COBRApy's gap-filling. Gap-filling needs a *universal* model — a pool of
candidate reactions to draw from (e.g. a BiGG universal model). Given a target
(grow, or produce a chosen metabolite), it returns the smallest set of reactions
from the universal model that makes the target feasible.
"""

from __future__ import annotations

from typing import List, Optional

import cobra
from cobra.flux_analysis import gapfill


class GapfillError(Exception):
    """Raised when gap-filling cannot be performed or finds no solution."""


def gapfill_for_growth(
    model: cobra.Model,
    universal: cobra.Model,
    *,
    lower_bound: float = 0.05,
    iterations: int = 1,
) -> List[List[str]]:
    """Find reaction sets from ``universal`` that let the model reach ``lower_bound`` growth."""
    try:
        solutions = gapfill(model, universal, lower_bound=lower_bound,
                            demand_reactions=False, iterations=iterations)
    except Exception as exc:  # noqa: BLE001
        raise GapfillError(f"Gap-filling failed:\n{exc}") from exc
    return [[r.id for r in sol] for sol in solutions]


def gapfill_for_metabolite(
    model: cobra.Model,
    universal: cobra.Model,
    metabolite_id: str,
    *,
    lower_bound: float = 0.05,
    iterations: int = 1,
) -> List[List[str]]:
    """Find reactions from ``universal`` that let the model *produce* ``metabolite_id``.

    A temporary demand reaction for the metabolite is set as the objective.
    """
    work = model.copy()
    if not work.metabolites.has_id(metabolite_id):
        raise GapfillError(f"No metabolite '{metabolite_id}' in the model.")
    met = work.metabolites.get_by_id(metabolite_id)
    demand = work.add_boundary(met, type="demand") if not _has_demand(work, met) else _get_demand(work, met)
    work.objective = demand
    try:
        solutions = gapfill(work, universal, lower_bound=lower_bound,
                            demand_reactions=False, iterations=iterations)
    except Exception as exc:  # noqa: BLE001
        raise GapfillError(
            f"No gap-filling solution found to produce {metabolite_id}.\n{exc}") from exc
    return [[r.id for r in sol] for sol in solutions]


def _has_demand(model: cobra.Model, met: cobra.Metabolite) -> bool:
    return any(r.id.startswith("DM_") and met in r.metabolites for r in model.reactions)


def _get_demand(model: cobra.Model, met: cobra.Metabolite):
    for r in model.reactions:
        if r.id.startswith("DM_") and met in r.metabolites:
            return r
    return None
