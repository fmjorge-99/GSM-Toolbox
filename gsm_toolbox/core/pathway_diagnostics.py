"""Why does a designed pathway carry no flux?

A route can be found, be sound chemistry, and still show zero flux — for two very
different reasons that users cannot tell apart from the number 0 alone:

* **Not incentivised.** The design-time prediction maximises production (objective =
  a demand on the product). A plain FBA maximises *biomass*, and making the product
  does not help the cell grow, so the optimiser routes nothing through the pathway.
  The route is fine; the objective is simply asking a different question. This is the
  common case behind "the tool predicted 0.21 but FBA gives 0".
* **Blocked.** Some step cannot carry flux in *any* steady state (FVA max = 0),
  usually because an intermediate is a dead end — only produced or only consumed. Then
  the route genuinely cannot run as it stands.

Flux Variability Analysis separates them: it asks "could this reaction EVER carry
flux?", independent of what the objective happens to reward.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import cobra

_TOL = 1e-9

# Verdicts
OK = "ok"                          # pathway carries flux at the current objective
NOT_INCENTIVISED = "not_incentivised"   # capable, but the objective doesn't reward it
BLOCKED = "blocked"                # a step cannot carry flux in any steady state
INFEASIBLE = "infeasible"          # the model itself cannot be solved


@dataclass
class FluxDiagnosis:
    """Why a pathway does (or does not) carry flux, and what to do about it."""

    target: str = ""
    reaction_ids: List[str] = field(default_factory=list)
    verdict: str = OK
    summary: str = ""
    recommendation: str = ""
    max_production: float = float("nan")      # with the objective on the product
    fluxes: Dict[str, float] = field(default_factory=dict)      # at current objective
    capacity: Dict[str, Tuple[float, float]] = field(default_factory=dict)  # FVA
    blocked_steps: List[str] = field(default_factory=list)
    bottleneck: str = ""              # first step that cannot carry flux
    last_carrying: str = ""           # last step before the bottleneck that can
    blocking_metabolites: List[str] = field(default_factory=list)
    objective: str = ""

    @property
    def carries_flux(self) -> bool:
        return self.verdict == OK


@dataclass
class NotFoundReason:
    """Why no route to a target could be found. "0 pathways" is never an answer."""

    target: str = ""
    reason: str = ""              # machine-readable: see the constants below
    summary: str = ""
    recommendation: str = ""
    missing: List[str] = field(default_factory=list)   # compounds nothing can make
    chain: List[str] = field(default_factory=list)     # target <- ... <- dead branch


NOT_IN_DB = "not_in_database"
NO_PRODUCER = "no_producing_reaction"
PRECURSOR_UNREACHABLE = "precursor_unreachable"
ALREADY_NATIVE = "already_native"
REACHABLE = "reachable"


def explain_not_found(host: cobra.Model, db: cobra.Model, target_id: str,
                      *, max_report: int = 6) -> NotFoundReason:
    """Explain why `target_id` has no route, naming the compound that breaks the chain.

    "0 pathways" hides three completely different situations: the compound may not be
    in the database at all; it may be there with no reaction that PRODUCES it (a
    database stub — the compound is catalogued, its biosynthesis is not); or it may be
    producible in principle but depend on a precursor the host cannot reach. Only the
    third is a search-limit problem, and only that one is helped by allowing more
    steps — which is why raising the step limit so often changes nothing.
    """
    from . import pathway_search as ps

    r = NotFoundReason(target=target_id)
    met = ps._resolve_target(host, db, target_id)
    if met is None:
        r.reason = NOT_IN_DB
        r.summary = (f"'{target_id}' is not in the loaded reaction database(s) at all.")
        r.recommendation = ("Load a database that contains this compound, merge several "
                            "databases, or fetch it online. Raising the step limit "
                            "cannot help — the compound simply is not there.")
        return r

    canon = ps.Canonicalizer([host, db])
    currency = ps._currency_keys(canon, [host, db])
    native = ps._native_keys(canon, host, True) | currency
    target_key = canon.key(met)
    if target_key in native:
        r.reason = ALREADY_NATIVE
        r.summary = f"Your model already makes '{target_id}' — no heterologous steps needed."
        return r

    # Which compounds can be produced at all, starting from the host's natives?
    producers: Dict[str, List[Tuple[str, set]]] = {}
    for rxn in db.reactions:
        if rxn.boundary:
            continue
        subs = {canon.key(m) for m, c in rxn.metabolites.items() if c < 0}
        prods = {canon.key(m) for m, c in rxn.metabolites.items() if c > 0}
        sub_need = {s for s in subs if s not in currency}
        prod_need = {p for p in prods if p not in currency}
        for p in prod_need - sub_need:
            producers.setdefault(p, []).append((rxn.id, sub_need))
        if rxn.reversibility or (rxn.lower_bound < 0 < rxn.upper_bound):
            for s in sub_need - prod_need:
                producers.setdefault(s, []).append((rxn.id, prod_need))

    if target_key not in producers:
        r.reason = NO_PRODUCER
        r.summary = (f"'{target_id}' is in the database, but NO reaction produces it. "
                     f"The compound is catalogued; its biosynthesis is not.")
        r.recommendation = ("This is a gap in the database, not a search limit — no "
                            "number of steps can help. Merge in another database, or "
                            "fetch the missing reactions online.")
        return r

    # Forward closure: everything the host can reach.
    reachable = set(native)
    changed = True
    while changed:
        changed = False
        for c, plist in producers.items():
            if c in reachable:
                continue
            for _rid, need in plist:
                if need <= reachable:
                    reachable.add(c)
                    changed = True
                    break
    if target_key in reachable:
        r.reason = REACHABLE
        r.summary = (f"'{target_id}' IS reachable — the search should find it. If it did "
                     f"not, the route may exceed the step limit.")
        r.recommendation = "Increase the maximum pathway length and search again."
        return r

    # Walk back from the target to the compounds that break the chain. A compound is a
    # blocker if it is unreachable AND none of its producing reactions can fire — i.e.
    # every producer still needs at least one OTHER unreachable non-currency compound.
    # This catches true orphans (no producer at all) AND circular pairs: cyclopentanol
    # ⇌ cyclopentanone each "produce" the other, but neither can be made from the host,
    # so cyclopentanone is the blocker to name even though it has a (circular) producer.
    seen: set = set()
    stack = [target_key]
    chain: List[str] = []
    while stack and len(r.missing) < max_report:
        c = stack.pop()
        if c in seen or c in reachable:
            continue
        seen.add(c)
        plist = producers.get(c)
        if not plist:
            r.missing.append(_pretty(db, canon, c))       # a true orphan
            continue
        # Can any producer fire from what's reachable, ignoring c itself? If a producer
        # needs only reachable compounds it would already be in `reachable` — so every
        # producer needs some unreachable substrate. Push those to keep walking; but if
        # every unreachable need is already `seen` (a closed circular/self-referential
        # cluster), c is the effective blocker.
        pushed = False
        for _rid, need in plist:
            for s in sorted(need):
                if s not in reachable and s != c and s not in seen:
                    stack.append(s)
                    pushed = True
        if not pushed:
            r.missing.append(_pretty(db, canon, c))       # unreachable circular cluster
        else:
            chain.append(_pretty(db, canon, c))

    r.reason = PRECURSOR_UNREACHABLE
    r.chain = chain[:max_report]
    if r.missing:
        names = ", ".join(r.missing[:4])
        r.summary = (f"A route to '{target_id}' exists on paper, but it depends on "
                     f"{names}, which NO reaction in this database produces. The chain "
                     f"cannot be closed at any step limit.")
        r.recommendation = (
            "This is a database gap. Merge another database (a compound catalogued in "
            "one is often fully synthesised in another), or fetch the missing reactions "
            "online. Raising the step limit will not help.")
    else:
        r.summary = (f"'{target_id}' cannot be reached from the compounds your model "
                     f"makes, using this database.")
        r.recommendation = ("Try merging additional databases, or allow more starting "
                            "metabolites.")
    return r


def _pretty(db: cobra.Model, canon, key: str) -> str:
    """A human-readable name for a canonical compound key."""
    for m in db.metabolites:
        if canon.key(m) == key:
            return f"{m.name or m.id} ({m.id})" if m.name else m.id
    return key


@dataclass
class BranchingAnalysis:
    """How much of the theoretical yield the surrounding network takes away.

    Following EA-MNE (J. Chem. Inf. Model. 2026, acs.jcim.5c02219): a designed route is
    not a line through an empty space. Every intermediate sits in the host's network,
    where other reactions compete for it. Those *branching* reactions divert carbon, so
    the yield a linear route promises is an upper bound the network may not deliver.
    """

    target: str = ""
    linear_yield: float = float("nan")     # max production, branches allowed to run
    network_yield: float = float("nan")    # max production with branches shut off
    branching_reactions: List[str] = field(default_factory=list)
    per_intermediate: Dict[str, List[str]] = field(default_factory=dict)
    yield_loss: float = float("nan")       # fraction of yield attributable to branching
    summary: str = ""

    @property
    def n_effective_branches(self) -> int:
        return len(self.branching_reactions)


def analyse_branching(model: cobra.Model, reaction_ids: List[str], *,
                      target_id: str) -> BranchingAnalysis:
    """Find the reactions competing with a designed route, and what they cost.

    `model` is the ENGINEERED model and is not modified. "Effective" branches are only
    those that can actually carry flux (FVA), which is the distinction EA-MNE draws:
    a competing reaction that is itself blocked costs nothing.
    """
    b = BranchingAnalysis(target=target_id)
    present = [r for r in reaction_ids if model.reactions.has_id(r)]
    if not present or not model.metabolites.has_id(target_id):
        return b
    route = set(present)

    def _max_production(m: cobra.Model) -> float:
        """Maximise production of the target.

        The sink must be resolved across compartments: `apply_pathway` typically adds a
        transport to the extracellular space and puts the exchange on the *_e* twin, so
        looking only at the cytosolic metabolite finds no boundary and the yield comes
        back null. Follow the transport, then fall back to a fresh demand (L11).
        """
        met = m.metabolites.get_by_id(target_id)
        sink = next((r for r in met.reactions if r.boundary), None)
        if sink is None:
            base = target_id.rsplit("_", 1)[0]
            for comp in ("e", "p"):
                twin = f"{base}_{comp}"
                if m.metabolites.has_id(twin):
                    sink = next((r for r in m.metabolites.get_by_id(twin).reactions
                                 if r.boundary), None)
                    if sink is not None:
                        break
        if sink is None:
            sink = m.add_boundary(met, type="demand")
        m.objective = sink
        v = m.slim_optimize()
        return float(v) if v is not None and v == v else float("nan")

    # Intermediates of the route (its non-currency compounds), and who else CONSUMES
    # them. Only consumers compete: a reaction that PRODUCES an intermediate is the
    # route's supply line, not a rival — including producers here would "prove" that
    # disabling the branches drops the yield to zero, which is just cutting the feed.
    from . import physiology
    from .pathway_design import _is_currency_met

    def _can_consume(r: cobra.Reaction, met: cobra.Metabolite) -> bool:
        c = r.metabolites[met]
        return (c < 0 and r.upper_bound > 1e-9) or (c > 0 and r.lower_bound < -1e-9)

    competitors: Dict[str, List[str]] = {}
    for rid in present:
        rxn = model.reactions.get_by_id(rid)
        for met, coeff in rxn.metabolites.items():
            # Only compounds the route CONSUMES can be diverted away from it.
            if coeff >= 0 or _is_currency_met(met) or met.id == target_id:
                continue
            # Biomass is excluded on purpose: it consumes nearly every precursor, so it
            # appears against every intermediate while being useless as a knock-down
            # target — the cell has to grow.
            others = [r.id for r in met.reactions
                      if r.id not in route and not r.boundary and _can_consume(r, met)
                      and not physiology.is_biomass_reaction(r)]
            if others:
                competitors.setdefault(met.id, [])
                for o in others:
                    if o not in competitors[met.id]:
                        competitors[met.id].append(o)

    all_others = sorted({o for lst in competitors.values() for o in lst})
    if not all_others:
        b.summary = "No reaction competes with this route for its intermediates."
        try:
            with model as m:
                b.linear_yield = b.network_yield = _max_production(m)
        except Exception:  # noqa: BLE001
            pass
        b.yield_loss = 0.0
        return b

    # Only branches that CAN carry flux divert anything.
    try:
        from cobra.flux_analysis import flux_variability_analysis
        fva = flux_variability_analysis(
            model, reaction_list=[model.reactions.get_by_id(r) for r in all_others],
            fraction_of_optimum=0.0, processes=1)   # spawn-safe (see producible_keys)
        effective = [r for r in all_others
                     if abs(float(fva.at[r, "minimum"])) > 1e-7
                     or abs(float(fva.at[r, "maximum"])) > 1e-7]
    except Exception:  # noqa: BLE001
        effective = all_others

    b.branching_reactions = effective
    b.per_intermediate = {m: [o for o in lst if o in set(effective)]
                          for m, lst in competitors.items()}
    b.per_intermediate = {m: lst for m, lst in b.per_intermediate.items() if lst}

    # "network" = the realistic figure, with the surrounding network free to compete.
    # "linear" = the ideal the route promises if nothing diverted its intermediates.
    # The gap between them is what branching costs.
    try:
        with model as m:
            b.network_yield = _max_production(m)
        with model as m:
            for r in effective:
                m.reactions.get_by_id(r).bounds = (0.0, 0.0)
            b.linear_yield = _max_production(m)
    except Exception:  # noqa: BLE001
        pass

    ly, ny = b.linear_yield, b.network_yield
    essential = (ly != ly or ly <= 1e-9) and (ny == ny and ny > 1e-9)
    if essential:
        # Disabling the branches broke the model: they are not optional side-drains but
        # reactions the cell needs (biomass precursors, cofactor recycling). Saying
        # "100% yield loss" here would be nonsense — the honest reading is that these
        # branches cannot simply be knocked out.
        b.yield_loss = float("nan")
        b.summary = (
            f"{len(effective)} reaction(s) compete for this route's intermediates, but "
            f"they cannot simply be removed — disabling them makes the model infeasible, "
            f"so they carry essential traffic (biomass precursors or cofactor recycling). "
            f"Maximum production with the network intact is {ny:.4g}. Down-regulation, "
            f"not knockout, is the lever here.")
        return b
    if ly == ly and ny == ny and ly > 1e-9:
        b.yield_loss = max(0.0, (ly - ny) / ly)
    b.summary = (
        f"{len(effective)} reaction(s) compete with this route for its intermediates. "
        f"Maximum production is {ny:.4g} with the network intact, and {ly:.4g} if those "
        f"branches were disabled.")
    if b.yield_loss == b.yield_loss and b.yield_loss > 0.01:
        b.summary += (f" About {b.yield_loss * 100:.0f}% of the achievable yield is lost "
                      f"to branching — these are the natural knockout / down-regulation "
                      f"targets.")
    elif b.yield_loss == b.yield_loss:
        b.summary += (" Branching costs almost nothing here: the optimiser can reach the "
                      "same yield without diverting carbon into them.")
    return b


def _dead_end_metabolites(model: cobra.Model, rxn: cobra.Reaction,
                          native_ids: Optional[set] = None) -> List[str]:
    """Metabolites of `rxn` that can only ever be produced, or only consumed — the
    usual reason a step is hard-blocked."""
    out = []
    for m in rxn.metabolites:
        if native_ids is not None and m.id in native_ids:
            continue
        makes = consumes = False
        for r in m.reactions:
            c = r.metabolites[m]
            fwd, rev = r.upper_bound > _TOL, r.lower_bound < -_TOL
            if (c > 0 and fwd) or (c < 0 and rev):
                makes = True
            if (c < 0 and fwd) or (c > 0 and rev):
                consumes = True
            if makes and consumes:
                break
        if not (makes and consumes):
            out.append(m.id)
    return out


def diagnose(model: cobra.Model, reaction_ids: List[str], *,
             target_id: str = "", native_ids: Optional[set] = None) -> FluxDiagnosis:
    """Diagnose the flux through `reaction_ids` in `model` (pathway order).

    `model` is the ENGINEERED model (host + pathway) and is not modified.
    """
    present = [r for r in reaction_ids if model.reactions.has_id(r)]
    d = FluxDiagnosis(target=target_id, reaction_ids=list(present))
    if not present:
        d.verdict = BLOCKED
        d.summary = "None of the pathway's reactions are present in the model."
        d.recommendation = "Re-apply the pathway to the model."
        return d

    try:
        d.objective = str(model.objective.expression)[:120]
    except Exception:  # noqa: BLE001
        d.objective = ""

    rxns = [model.reactions.get_by_id(r) for r in present]

    # 1) Can each step carry flux at all? FVA answers this independently of the
    #    objective, which is exactly the distinction the user cannot see from "0".
    try:
        from cobra.flux_analysis import flux_variability_analysis
        fva = flux_variability_analysis(model, reaction_list=rxns, fraction_of_optimum=0.0)
        for rid in present:
            d.capacity[rid] = (float(fva.at[rid, "minimum"]), float(fva.at[rid, "maximum"]))
    except Exception as exc:  # noqa: BLE001
        d.verdict = INFEASIBLE
        d.summary = f"The model could not be analysed: {exc}"
        d.recommendation = "Check the model solves (run FBA) before diagnosing the pathway."
        return d

    d.blocked_steps = [rid for rid, (lo, hi) in d.capacity.items()
                       if abs(lo) <= 1e-7 and abs(hi) <= 1e-7]

    # 2) What does the CURRENT objective actually do with the pathway?
    try:
        sol = model.optimize()
        if sol.status == "optimal":
            d.fluxes = {rid: float(sol.fluxes[rid]) for rid in present}
    except Exception:  # noqa: BLE001
        pass

    # 3) Maximum production, for the "what could it do" figure.
    if target_id and model.metabolites.has_id(target_id):
        try:
            with model as m:
                met = m.metabolites.get_by_id(target_id)
                sink = None
                for r in met.reactions:
                    if r.boundary:
                        sink = r
                        break
                if sink is None:
                    sink = m.add_boundary(met, type="demand")
                m.objective = sink
                val = m.slim_optimize()
                d.max_production = float(val) if val is not None and val == val else float("nan")
        except Exception:  # noqa: BLE001
            pass

    # --- verdict -------------------------------------------------------------
    if d.blocked_steps:
        d.verdict = BLOCKED
        # The bottleneck is the FIRST step (in pathway order) that cannot carry flux;
        # everything before it is fine, which is what the user needs to know.
        order = {rid: i for i, rid in enumerate(present)}
        d.bottleneck = min(d.blocked_steps, key=lambda r: order[r])
        idx = order[d.bottleneck]
        carrying = [r for r in present[:idx] if r not in d.blocked_steps]
        d.last_carrying = carrying[-1] if carrying else ""
        # The dead end that blocks a step is often DOWNSTREAM of it: an intermediate
        # further along has no consumer, which back-propagates and pins the earlier
        # steps to zero too. So look across every blocked step, not just the first.
        seen = set()
        for rid in sorted(d.blocked_steps, key=lambda r: order[r]):
            for mid in _dead_end_metabolites(model, model.reactions.get_by_id(rid),
                                             native_ids):
                if mid not in seen:
                    seen.add(mid)
                    d.blocking_metabolites.append(mid)
        where = (f" The last step that can carry flux is {d.last_carrying}; "
                 f"the route stops at {d.bottleneck}."
                 if d.last_carrying else f" The route is blocked at its first step, "
                                         f"{d.bottleneck}.")
        d.summary = (f"{len(d.blocked_steps)} step(s) cannot carry flux in any steady "
                     f"state, so this pathway cannot run as it stands.{where}")
        if d.blocking_metabolites:
            names = ", ".join(d.blocking_metabolites[:4])
            d.recommendation = (
                f"{d.bottleneck} is blocked because {names} can only be produced or only "
                f"consumed (a dead end). Add a reaction that consumes/produces it, give "
                f"it an exchange if it should leave the cell, or pick an alternative "
                f"route that avoids this step.")
        else:
            d.recommendation = (
                f"{d.bottleneck} cannot carry flux — check its bounds and that all of "
                f"its substrates can be made by the host.")
        return d

    if d.fluxes and all(abs(v) <= 1e-7 for v in d.fluxes.values()):
        d.verdict = NOT_INCENTIVISED
        cap = max((hi for _lo, hi in d.capacity.values()), default=0.0)
        d.summary = (
            "Every step CAN carry flux (up to "
            f"{cap:.4g}), but at the current objective the pathway carries none. "
            "This is expected: FBA maximises the objective, and making the product "
            "does not help it — so the optimiser sends no carbon down the route. "
            "The pathway is not broken.")
        d.recommendation = (
            "To see production, either set the objective to the product's exchange/demand "
            "reaction (Maximize Product), or keep biomass as the objective and force a "
            "minimum production by raising the lower bound of the product exchange. "
            "A production envelope shows the achievable growth/production trade-off.")
        return d

    d.verdict = OK
    carried = max((abs(v) for v in d.fluxes.values()), default=0.0)
    d.summary = f"The pathway carries flux at the current objective (up to {carried:.4g})."
    d.recommendation = ""
    return d
