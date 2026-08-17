"""Which metabolites the host *actually makes* while growing — not merely which it lists.

A genome-scale model contains many metabolites that are present in the stoichiometry but
carry no flux under the growth conditions being simulated. The pathway search treats every
native compound as a valid place for a route to begin, so a design can be reported as
"starting from a native precursor" when the host does not, in fact, produce that precursor
at all.

**Lactaldehyde in Synechocystis PCC 6803 is the worked example.** It appears in the model,
so a 1,2-propanediol route can end its search there and look complete — but wild-type
Synechocystis has no route to lactaldehyde. A real strain needs a heterologous
methylglyoxal synthase (*mgsA*) or glycerol dehydrogenase (*gldA*) to supply it. A design
that quietly assumes lactaldehyde is available is not buildable, and nothing in the output
says so.

This module answers the question that distinguishes the two cases: *under the growth
conditions in the model right now, which metabolites carry flux?* Restricting the search's
starting pool to those compounds is all it takes to force the search to bridge the gap
itself — it must then find the heterologous reactions (the *mgsA* step) that connect a
genuinely available precursor to the target, because it is no longer allowed to start from
one that is idle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Set

import cobra

#: Below this a flux is numerical noise rather than a real carrier.
FLUX_TOL = 1e-6


@dataclass
class FluxContext:
    """The metabolites a host produces and consumes at a growth solution."""

    objective_value: float = 0.0
    status: str = ""
    produced: Set[str] = field(default_factory=set)
    consumed: Set[str] = field(default_factory=set)
    active_reactions: Set[str] = field(default_factory=set)
    #: Native metabolites present in the model but carrying no flux at all.
    idle: Set[str] = field(default_factory=set)

    @property
    def usable(self) -> bool:
        """True when the solution is meaningful enough to restrict a search by."""
        return self.status == "optimal" and bool(self.produced)

    def summary(self) -> str:
        if self.status != "optimal":
            return (f"The model could not be solved ({self.status or 'no solution'}), so "
                    "flux-carrying precursors cannot be identified.")
        if self.objective_value <= FLUX_TOL:
            return ("The objective is zero at these conditions — nothing is growing, so "
                    "no metabolite carries flux. Check the medium and the objective.")
        return (f"{len(self.produced)} metabolite(s) are actively produced at the current "
                f"growth solution (objective {self.objective_value:.4g}); "
                f"{len(self.idle)} native metabolite(s) carry no flux and are therefore "
                "not available as starting points.")


def growth_flux_context(model: cobra.Model, *, parsimonious: bool = True,
                        tol: float = FLUX_TOL) -> FluxContext:
    """Solve the model as configured and report which metabolites carry flux.

    ``parsimonious`` uses pFBA, which removes the arbitrary internal loops a plain FBA
    solution can contain — those loops would otherwise mark metabolites as "produced"
    when nothing real is making them. It falls back to plain FBA if pFBA is unavailable
    or infeasible, because a slightly loose answer is far more useful than none.
    """
    ctx = FluxContext()
    try:
        solution = None
        if parsimonious:
            try:
                from cobra.flux_analysis import pfba
                solution = pfba(model)
            except Exception:  # noqa: BLE001 - loopless solve unavailable/infeasible
                solution = None
        if solution is None:
            solution = model.optimize()
        ctx.status = str(solution.status)
        ctx.objective_value = float(solution.objective_value or 0.0)
    except Exception as exc:  # noqa: BLE001 - an unsolvable model is a normal outcome
        ctx.status = f"error: {exc}"
        return ctx
    if ctx.status != "optimal":
        return ctx

    fluxes = solution.fluxes
    for rxn in model.reactions:
        v = float(fluxes.get(rxn.id, 0.0) or 0.0)
        if abs(v) < tol:
            continue
        ctx.active_reactions.add(rxn.id)
        for met, coeff in rxn.metabolites.items():
            # Net direction at this solution, not the direction as written.
            net = coeff * v
            if net > tol:
                ctx.produced.add(met.id)
            elif net < -tol:
                ctx.consumed.add(met.id)
    touched = ctx.produced | ctx.consumed
    ctx.idle = {m.id for m in model.metabolites} - touched
    return ctx


def flux_carrying_starts(model: cobra.Model, *, requested: Optional[list] = None,
                         parsimonious: bool = True) -> tuple:
    """``(start_ids, context)`` — the precursors a route may legitimately begin from.

    When ``requested`` is given (the user picked specific starting metabolites), it is
    *intersected* with the flux-carrying set rather than replaced: an explicit choice
    still narrows the search, but a compound the host is not making cannot sneak back in.
    """
    ctx = growth_flux_context(model, parsimonious=parsimonious)
    if not ctx.usable:
        # Never silently return an empty pool — that would make every search fail. The
        # caller reports the problem and falls back to the unrestricted behaviour.
        return (list(requested or []), ctx)
    starts = set(ctx.produced)
    if requested:
        chosen = set(requested)
        overlap = chosen & starts
        # If the user's picks are all idle, honour the picks and let the caller warn:
        # overriding them silently would be worse than an explicit, explained conflict.
        starts = overlap or chosen
    return (sorted(starts), ctx)


def idle_native_precursors(model: cobra.Model, metabolite_ids, *,
                           context: Optional[FluxContext] = None) -> Dict[str, str]:
    """Of ``metabolite_ids``, which are present in the host but carry no flux.

    Used to explain a route: "this design begins at lactaldehyde, which your host does
    not currently make".
    """
    ctx = context or growth_flux_context(model)
    out: Dict[str, str] = {}
    if ctx.status != "optimal":
        return out
    for mid in metabolite_ids:
        if not model.metabolites.has_id(mid):
            continue
        if mid in ctx.produced:
            continue
        met = model.metabolites.get_by_id(mid)
        out[mid] = met.name or mid
    return out
