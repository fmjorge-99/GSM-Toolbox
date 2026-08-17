"""Mutant-phenotype simulation: MOMA and ROOM.

After a gene/reaction knockout a cell does not immediately re-optimize growth like
FBA assumes; it tends to stay close to its wild-type flux state. MOMA and ROOM
capture this and usually predict knockout phenotypes better than plain FBA:

* **MOMA** (Minimization Of Metabolic Adjustment) — the mutant flux distribution
  that is closest to the wild-type one [Segre et al. 2002].
* **ROOM** (Regulatory On/Off Minimization) — minimizes the *number* of reactions
  whose flux changes significantly from wild-type [Shlomi et al. 2005].
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import cobra
import pandas as pd
from cobra.flux_analysis import moma, pfba, room

from ..editing import guess_biomass_reaction


class MutantError(Exception):
    """Raised when a mutant simulation cannot be performed."""


@dataclass
class MutantResult:
    method: str
    status: str
    growth: float
    fluxes: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    wt_growth: float = float("nan")
    n_changed: int = 0

    def flux_table(self) -> pd.DataFrame:
        df = self.fluxes.rename("mutant_flux").to_frame()
        df.index.name = "reaction"
        return df.reset_index()


def run_mutant(
    model: cobra.Model,
    knockouts: List[str],
    *,
    method: str = "moma",
    biomass_id: Optional[str] = None,
) -> MutantResult:
    """Predict a knockout mutant's flux distribution with MOMA or ROOM.

    ``method`` is 'moma' (linear MOMA) or 'room'. ``knockouts`` are reaction ids to
    delete. The wild-type reference is the parsimonious FBA solution.
    """
    if not knockouts:
        raise MutantError(
            "Select one or more reactions to knock out (in the Explorer) before running "
            "a mutant simulation.")
    biomass_id = biomass_id or guess_biomass_reaction(model)

    try:
        wt = pfba(model)
    except Exception as exc:  # noqa: BLE001
        raise MutantError(f"Wild-type reference optimization failed:\n{exc}") from exc
    wt_growth = float(wt.fluxes.get(biomass_id, float("nan"))) if biomass_id else float("nan")

    with model:
        applied = []
        for rid in knockouts:
            if model.reactions.has_id(rid):
                model.reactions.get_by_id(rid).knock_out()
                applied.append(rid)
        if not applied:
            raise MutantError("None of the selected reactions exist in this model.")
        try:
            if method == "room":
                sol = room(model, solution=wt, linear=True)
            else:
                sol = moma(model, solution=wt, linear=True)
        except Exception as exc:  # noqa: BLE001
            raise MutantError(
                f"{method.upper()} failed:\n{exc}\n(quadratic MOMA needs a QP solver such "
                "as CPLEX/Gurobi; linear MOMA and ROOM work with the bundled solvers.)") from exc

    growth = float(sol.fluxes.get(biomass_id, float("nan"))) if biomass_id else float("nan")
    # Reactions whose flux changed appreciably from wild-type.
    diff = (sol.fluxes - wt.fluxes).abs()
    n_changed = int((diff > 1e-6).sum())
    return MutantResult(
        method=method.upper(), status=sol.status, growth=growth, fluxes=sol.fluxes,
        wt_growth=wt_growth, n_changed=n_changed)
