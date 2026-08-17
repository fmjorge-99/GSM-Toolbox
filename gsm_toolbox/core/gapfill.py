"""Find what actually blocks a designed pathway, and close the gap automatically.

Repeatedly across real targets — violacein, THC, psilocybin, n-butanol — the toolbox found
a chemically sound route that carried no flux, and the honest diagnosis was always the
same shape: *one metabolite in the route cannot be supplied by the host*. Violacein's
first oxidation makes an intermediate nothing consumes; n-butanol's β-oxidation step needs
butyryl-CoA, which *Synechocystis* does not make.

Telling the user that is useful. Doing something about it is far more useful. This module:

1. **Names the blocking metabolite**, not just the blocked reaction. A topological
   dead-end test is not enough — butyryl-CoA *has* producers in the merged database, they
   simply cannot carry flux — so the test here is the operational one: *can this compound
   be produced at all?* (an FVA-style max-production check).
2. **Fetches the missing chemistry around it** from the open sources (KEGG, Rhea) and
   re-searches, returning a complete route in one step rather than asking the user to
   assemble it by hand.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import cobra

_TOL = 1e-7


@dataclass
class Blocker:
    """A metabolite the host cannot supply, and where it hurts."""

    metabolite_id: str
    name: str = ""
    reaction_id: str = ""          # the step that needs it
    reason: str = ""               # human-readable cause

    @property
    def label(self) -> str:
        return f"{self.name} ({self.metabolite_id})" if self.name else self.metabolite_id


@dataclass
class BlockageReport:
    """Why a route carries no flux, in terms the user can act on."""

    blocked: bool = False
    bottleneck: str = ""                        # first step that cannot carry flux
    blockers: List[Blocker] = field(default_factory=list)
    recommendation: str = ""

    @property
    def blocker_labels(self) -> List[str]:
        return [b.label for b in self.blockers]

    def sentence(self) -> str:
        if not self.blocked:
            return "This route can carry flux."
        if self.blockers:
            names = ", ".join(b.label for b in self.blockers[:3])
            return (f"Blocked at {self.bottleneck}: the host cannot supply {names}. "
                    f"{self.recommendation}")
        return (f"Blocked at {self.bottleneck}, but no single missing precursor could be "
                f"identified. {self.recommendation}")


def _can_produce(model: cobra.Model, met_id: str) -> bool:
    """Can the model make this metabolite at all? (max production > 0)"""
    if not model.metabolites.has_id(met_id):
        return False
    try:
        with model as m:
            met = m.metabolites.get_by_id(met_id)
            sink = next((r for r in met.reactions if r.boundary), None)
            if sink is None:
                sink = m.add_boundary(met, type="demand")
            m.objective = sink
            v = m.slim_optimize()
            return bool(v is not None and v == v and abs(v) > 1e-9)
    except Exception:  # noqa: BLE001
        return False


def _can_consume(model: cobra.Model, met_id: str) -> bool:
    """Can the model dispose of this metabolite? (a sink can carry flux)"""
    if not model.metabolites.has_id(met_id):
        return False
    try:
        with model as m:
            met = m.metabolites.get_by_id(met_id)
            sink = m.add_boundary(met, type="sink")
            m.objective = sink
            m.objective_direction = "min"      # a sink consuming runs negative
            v = m.slim_optimize()
            return bool(v is not None and v == v and abs(v) > 1e-9)
    except Exception:  # noqa: BLE001
        return False


_CURRENCY_HINTS = ("h2o", "water", "h_c", "h_p", "h_e", "proton", "co2", "o2", "nad",
                   "nadp", "fad", "fmn", "atp", "adp", "amp", "pi_", "ppi", "coa_c",
                   "nh4", "nh3", "h2o2", "photon")


def _is_currency(met) -> bool:
    mid = (getattr(met, "id", "") or "").lower()
    name = (getattr(met, "name", "") or "").lower()
    if mid.rsplit("_", 1)[0] in ("h", "h2o", "co2", "o2", "pi", "ppi", "nh4", "h2o2"):
        return True
    return any(h in mid for h in _CURRENCY_HINTS) or name in (
        "h+", "h2o", "water", "co2", "o2", "oxygen", "phosphate", "diphosphate")


def find_blockers(model: cobra.Model, reaction_ids: List[str], *,
                  native_ids: Optional[set] = None) -> BlockageReport:
    """Identify the bottleneck step AND the metabolites that make it impossible.

    The key difference from a topological dead-end test: a compound can have producers on
    paper yet still be unmakeable, which is exactly the butyryl-CoA case. Every substrate
    of the blocked step is therefore tested by asking the model to maximise its production.
    """
    rep = BlockageReport()
    present = [r for r in reaction_ids if model.reactions.has_id(r)]
    if not present:
        return rep

    # Which steps can carry flux at all?
    blocked_steps: List[str] = []
    try:
        from cobra.flux_analysis import flux_variability_analysis
        fva = flux_variability_analysis(
            model, reaction_list=[model.reactions.get_by_id(r) for r in present],
            fraction_of_optimum=0.0, processes=1)
        for rid in present:
            lo, hi = float(fva.at[rid, "minimum"]), float(fva.at[rid, "maximum"])
            if abs(lo) <= _TOL and abs(hi) <= _TOL:
                blocked_steps.append(rid)
    except Exception:  # noqa: BLE001
        return rep
    if not blocked_steps:
        return rep

    rep.blocked = True
    order = {rid: i for i, rid in enumerate(present)}
    rep.bottleneck = min(blocked_steps, key=lambda r: order[r])

    # Which metabolite actually makes the route impossible? Two failure modes, and the
    # culprit is often NOT in the first blocked step: a compound with no consumer several
    # steps downstream backs the whole chain up, pinning the earlier steps to zero. So
    # examine every blocked step, in pathway order.
    #
    #   * a SUBSTRATE that cannot be produced  → the step starves;
    #   * a PRODUCT that cannot be consumed    → the step backs up.
    #
    # Note a reversible reaction can often make its own substrate from the other side, so
    # "cannot be produced" is a genuine test, not a formality.
    seen: set = set()
    for rid in sorted(blocked_steps, key=lambda r: order[r]):
        rxn = model.reactions.get_by_id(rid)
        for met, coeff in rxn.metabolites.items():
            if met.id in seen:
                continue
            if _is_currency(met):
                continue                  # water, protons, NAD… never the real blocker
            if coeff < 0 and not _can_produce(model, met.id):
                seen.add(met.id)
                rep.blockers.append(Blocker(
                    metabolite_id=met.id, name=getattr(met, "name", "") or "",
                    reaction_id=rid,
                    reason="nothing in the model can produce it, so this step starves"))
            elif coeff > 0 and not _can_consume(model, met.id):
                seen.add(met.id)
                rep.blockers.append(Blocker(
                    metabolite_id=met.id, name=getattr(met, "name", "") or "",
                    reaction_id=rid,
                    reason="nothing in the model consumes it, so this step backs up"))

    if rep.blockers:
        first = rep.blockers[0]
        rep.recommendation = (
            f"Add chemistry that produces {first.label} — either fetch reactions around "
            f"it from KEGG/Rhea, or design a separate pathway to it and add that first.")
    else:
        rep.recommendation = (
            f"Check the bounds of {rep.bottleneck} and whether its by-products can be "
            "disposed of; nothing it consumes is individually unmakeable.")
    return rep


# --------------------------------------------------------------------------------------
# Automatic gap closing
# --------------------------------------------------------------------------------------
@dataclass
class GapFillResult:
    """Outcome of an automatic gap-closing attempt."""

    success: bool = False
    fetched: List[dict] = field(default_factory=list)   # {source, term, n_reactions}
    reactions_added: int = 0
    reaction_ids: List[str] = field(default_factory=list)   # the completed route
    target_id: str = ""
    max_flux: Optional[float] = None
    message: str = ""
    merged_model: Optional[cobra.Model] = None
    still_blocked: List[str] = field(default_factory=list)


def autofill(host: cobra.Model, db: cobra.Model, target_id: str,
             blockers: List[Blocker], *, sources: Tuple[str, ...] = ("KEGG", "Rhea"),
             expand_steps: int = 1, max_reactions: int = 200,
             progress: Optional[Callable[[str], None]] = None) -> GapFillResult:
    """Fetch chemistry around the blocking metabolites and re-search for a complete route.

    Returns a :class:`GapFillResult` whose ``reaction_ids`` and ``merged_model`` are a
    route the caller can apply directly — the point being that the user gets a finished
    pathway rather than a shopping list.
    """
    from . import databases, pathway_design as pdz

    res = GapFillResult(target_id=target_id)
    if not blockers:
        res.message = "No blocking metabolite was identified, so there is nothing to fetch."
        return res

    extra: List[cobra.Model] = []
    for b in blockers[:3]:                       # a few blockers at most; keep it quick
        term = (b.name or b.metabolite_id).strip()
        # Strip a compartment suffix from a bare id so the online search has a chance.
        if not b.name and "_" in term:
            term = term.rsplit("_", 1)[0]
        for src in sources:
            if progress:
                progress(f"Fetching {src} reactions around “{term}”…")
            try:
                fn = (databases.build_kegg_pathway_db if src == "KEGG"
                      else databases.build_rhea_pathway_db)
                model, label, _p = fn(term, expand_steps=expand_steps,
                                      max_reactions=max_reactions)
                res.fetched.append({"source": src, "term": term, "label": label,
                                    "n_reactions": len(model.reactions)})
                extra.append(model)
            except Exception as exc:  # noqa: BLE001 — a missing source is not fatal
                res.fetched.append({"source": src, "term": term,
                                    "error": str(exc)[:160]})

    if not extra:
        res.message = ("Neither KEGG nor Rhea has chemistry around "
                       f"{blockers[0].label}. This gap cannot be closed from the open "
                       "sources; it needs a curated reaction or an operon transplant.")
        return res

    if progress:
        progress("Merging the fetched chemistry…")
    merged = db.copy()
    added = 0
    for m in extra:
        for r in m.reactions:
            if merged.reactions.has_id(r.id):
                continue
            try:
                merged.add_reactions([r.copy()])
                added += 1
            except Exception:  # noqa: BLE001
                pass
    res.reactions_added = added
    res.merged_model = merged

    if progress:
        progress(f"Re-searching with {added} new reactions…")
    try:
        found = pdz.find_pathways(host, target_id, merged, max_steps=16, n_alternatives=3)
        routes = [r for r in found if r.reaction_ids]
    except Exception as exc:  # noqa: BLE001
        res.message = f"The re-search failed: {exc}"
        return res

    if not routes:
        res.message = (f"Added {added} reactions around {blockers[0].label}, but still no "
                       "route to the target. The gap is deeper than one compound.")
        return res

    # Prefer a route that actually carries flux.
    best = max(routes, key=lambda r: (r.production_flux == r.production_flux
                                      and r.production_flux or 0.0))
    res.reaction_ids = list(best.reaction_ids)
    f = best.production_flux
    res.max_flux = None if f != f else float(f)
    res.success = True
    if res.max_flux and res.max_flux > 1e-9:
        res.message = (f"Fetched {added} reactions around {blockers[0].label} and found a "
                       f"complete {len(res.reaction_ids)}-step route carrying "
                       f"{res.max_flux:.4g} mmol gDW⁻¹ h⁻¹.")
    else:
        res.message = (f"Fetched {added} reactions around {blockers[0].label} and found a "
                       f"complete {len(res.reaction_ids)}-step route, but it still carries "
                       "no flux — run the flux diagnosis again to see what now blocks it.")
    return res
