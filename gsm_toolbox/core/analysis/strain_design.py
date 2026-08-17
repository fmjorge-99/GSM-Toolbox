"""Strain design for metabolic engineering (Phase 4).

Two complementary families:

* **Bi-level knockout design** (OptKnock / RobustKnock) via the StrainDesign
  package — finds reaction knockout sets that couple a target product to growth.
* **FSEOF** (Flux Scanning based on Enforced Objective Flux) — implemented here —
  scans increasing enforced product flux and classifies reactions as
  over-expression or knock-down/deletion targets based on their flux trend.

All functions operate on a model copy passed in by the caller and return tidy
result objects for the GUI.
"""

from __future__ import annotations

import contextlib
import logging
import warnings
from dataclasses import dataclass, field
from typing import List, Optional

import cobra
import numpy as np
import pandas as pd

from .. import physiology
from ..editing import guess_biomass_reaction


class StrainDesignError(Exception):
    """Raised when a strain-design computation cannot run or finds nothing."""


@dataclass
class StrainDesignResult:
    method: str
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    note: str = ""


# --- Bi-level knockout design (StrainDesign) --------------------------------
def _solver_available() -> str:
    """Report the best available MILP solver (affects feasibility of OptKnock)."""
    import straindesign as sd

    avail = sd.avail_solvers
    for pref in ("gurobi", "cplex", "scip", "glpk"):
        if pref in avail:
            return pref
    return "glpk"


def _is_transport(rxn: cobra.Reaction) -> bool:
    """Heuristic: a reaction is a transporter if a metabolite appears in >1 compartment."""
    bases = {}
    for met in rxn.metabolites:
        base = met.id.rsplit("_", 1)[0] if "_" in met.id else met.id
        bases.setdefault(base, set()).add(met.compartment)
    return any(len(comps) > 1 for comps in bases.values())


EXACT_METHODS = {"optknock", "robustknock", "optcouple"}

# GLPK is technically a MILP solver but is slow/unstable on genome-scale bi-level
# problems (endless "infeasible solution … a subset seems valid" churn, near-hangs).
# These are the solvers we consider *capable* of OptKnock-class problems.
_CAPABLE_MILP = ("gurobi", "cplex", "scip")


@contextlib.contextmanager
def _quiet_solver():
    """Aggregate/suppress the repetitive solver churn (Issue 8).

    StrainDesign/GLPK emit a stream of "Solver first found the infeasible
    solution … a subset seems valid" messages that flood the log. Silence the
    relevant loggers and Python warnings for the duration of the solve so the log
    stays readable; levels are restored afterward.
    """
    names = ["straindesign", "cobra", "optlang", "pyscipopt", "swiglpk"]
    loggers = [logging.getLogger(n) for n in names]
    prev = [lg.level for lg in loggers]
    for lg in loggers:
        lg.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            yield
    finally:
        for lg, lvl in zip(loggers, prev):
            lg.setLevel(lvl)


def _capable_milp_solver() -> Optional[str]:
    """Return the best genome-scale-capable MILP solver available, or None (Issue 8)."""
    try:
        import straindesign as sd
        avail = set(sd.avail_solvers)
    except Exception:  # noqa: BLE001
        return None
    for pref in _CAPABLE_MILP:
        if pref in avail:
            return pref
    return None


def available_solvers() -> List[str]:
    """Solver names StrainDesign can use right now (e.g. ['scip', 'glpk'])."""
    import straindesign as sd

    return sorted(sd.avail_solvers)


def secrete_metabolite(model: cobra.Model, met_id: str, *, gas: bool = False) -> str:
    """Ensure a realistic secretion route for a product metabolite (Issue 10).

    Rather than only draining the cytosolic species with a ``DM_`` demand, add
    (if missing) a transport step to the extracellular compartment and an
    exchange reaction — so production envelopes / titers reflect real secretion.
    For volatile products (``gas=True``) the exchange represents a gas-phase sink.
    Returns the id of the exchange (or demand) reaction to use as the product.
    """
    if not model.metabolites.has_id(met_id):
        raise StrainDesignError(f"No metabolite '{met_id}' in the model.")
    met = model.metabolites.get_by_id(met_id)
    comps = set(getattr(model, "compartments", {}) or {}) or {
        m.compartment for m in model.metabolites if m.compartment}
    # No extracellular compartment: fall back to a demand drain (best we can do).
    ext = "e" if "e" in comps else next((c for c in comps if c.startswith("e")), None)
    base = met.id.rsplit("_", 1)[0] if "_" in met.id else met.id
    if ext is None or met.compartment == ext:
        # Already extracellular (or no 'e'): just add an exchange/demand.
        ex_id = f"EX_{met.id}"
        if not model.reactions.has_id(ex_id):
            model.add_boundary(met, type="exchange")
        return ex_id

    ext_id = f"{base}_{ext}"
    if model.metabolites.has_id(ext_id):
        ext_met = model.metabolites.get_by_id(ext_id)
    else:
        ext_met = cobra.Metabolite(ext_id, name=met.name, formula=met.formula,
                                   charge=met.charge, compartment=ext)
        model.add_metabolites([ext_met])
    # Transport met_c -> met_e
    t_id = f"{base}t" + ("_gas" if gas else "")
    if not model.reactions.has_id(t_id):
        t = cobra.Reaction(t_id, name=f"{met.name or base} transport", lower_bound=0,
                           upper_bound=1000)
        model.add_reactions([t])
        t.add_metabolites({met: -1, ext_met: 1})
    # Exchange / gas-phase sink
    ex_id = f"EX_{ext_id}"
    if not model.reactions.has_id(ex_id):
        model.add_boundary(ext_met, type="exchange")
    return ex_id


# How a product is modelled as leaving (or staying in) the cell. In a stoichiometric
# model active secretion and passive diffusion are represented identically (a
# membrane transport step + an extracellular exchange); a volatile product uses the
# same route flagged as a gas-phase escape; an intracellularly-accumulated product
# gets a cytosolic demand sink (no membrane step); "none" adds nothing.
PRODUCT_ROUTES = ("secrete", "diffusion", "gas", "intracellular", "none")


def ensure_product_exchange(model: cobra.Model, met_id: str, *, route: str = "secrete",
                            gas: bool = False):
    """Guarantee the product has exactly one route out of (or sink within) the model,
    the way an expert would set it up before measuring a titre. If a boundary reaction
    already drains this compound it is reused (never a duplicate). ``route`` selects the
    mechanism (see :data:`PRODUCT_ROUTES`). Returns ``(reaction_id_or_None, created,
    route_kind)``."""
    if gas and route == "secrete":
        route = "gas"
    if route not in PRODUCT_ROUTES:
        route = "secrete"
    if route == "none":
        return None, False, "none"
    if not model.metabolites.has_id(met_id):
        raise StrainDesignError(f"No metabolite '{met_id}' in the model.")
    met = model.metabolites.get_by_id(met_id)
    base = met.id.rsplit("_", 1)[0] if "_" in met.id else met.id

    def _base(mid):
        return mid.rsplit("_", 1)[0] if "_" in mid else mid

    if route == "intracellular":
        # Product stays inside the cell: a cytosolic demand sink lets it be produced
        # at steady state (models storage/accumulation), with no membrane transport.
        for rxn in met.reactions:
            if rxn.boundary:
                return rxn.id, False, "intracellular"
        return model.add_boundary(met, type="demand").id, True, "intracellular"

    # secrete / diffusion / gas: reuse-or-create a transport + extracellular exchange.
    candidates = {met}
    for m in model.metabolites:
        if _base(m.id) == base and (m.compartment or "").startswith("e"):
            candidates.add(m)
    for m in candidates:
        for rxn in m.reactions:
            if rxn.boundary:
                return rxn.id, False, route
    return secrete_metabolite(model, met_id, gas=(route == "gas")), True, route


def _knockable_reactions(model, product_id, biomass_id, *,
                         exclude_exchanges=True, exclude_transport=False,
                         exclude_blocked=True) -> set:
    """Reaction ids that may be knocked out, after sensible exclusions."""
    excluded = {biomass_id, product_id}
    if model.reactions.has_id("ATPM"):
        excluded.add("ATPM")
    if exclude_exchanges:
        excluded |= {r.id for r in model.reactions if r.boundary}
    if exclude_transport:
        excluded |= {r.id for r in model.reactions if _is_transport(r)}
    if exclude_blocked:
        try:
            from cobra.flux_analysis import find_blocked_reactions
            excluded |= set(find_blocked_reactions(model))
        except Exception:  # noqa: BLE001
            pass
    return {r.id for r in model.reactions if r.id not in excluded}


def _wt_growth(model, biomass_id) -> float:
    with model:
        model.objective = biomass_id
        g = model.slim_optimize()
    return float(g) if g is not None and g == g else 0.0


def run_knockout_design(
    model: cobra.Model,
    product_id: str,
    *,
    method: str = "optknock",
    solver: str = "auto",
    **kwargs,
) -> StrainDesignResult:
    """Dispatch a knockout strain-design run to the exact (MILP) or heuristic engine.

    Issue 8: exact (OptKnock-class) methods need a capable MILP solver. If the
    user asked for 'auto' and only GLPK is present, fall back to the built-in
    evolutionary heuristic (which is far more robust here) rather than churning
    on GLPK — and say so in the result note.
    """
    if method == "heuristic":
        return heuristic_knockout_design(model, product_id, **kwargs)

    if method in EXACT_METHODS and (not solver or solver == "auto") \
            and _capable_milp_solver() is None:
        kwargs.pop("solver", None)
        kwargs.pop("approach", None)   # MILP-only; heuristic ignores extras anyway
        result = heuristic_knockout_design(model, product_id, **kwargs)
        result.note = ("No capable MILP solver (SCIP/Gurobi/CPLEX) was found, so the "
                       "built-in evolutionary heuristic was used instead of "
                       f"{method.upper()} (GLPK is unreliable on genome-scale bi-level "
                       "problems). " + result.note)
        return result
    try:
        return compute_knockout_design(model, product_id, method=method, solver=solver, **kwargs)
    except StrainDesignError:
        raise
    except Exception as exc:  # noqa: BLE001
        # The exact solver failed unexpectedly (e.g. RobustKnock's known SCIP quirk).
        # Fall back to the robust heuristic rather than leaving the user stuck.
        for k in ("solver", "approach"):
            kwargs.pop(k, None)
        result = heuristic_knockout_design(model, product_id, **kwargs)
        result.note = (f"{method.upper()} could not be solved exactly ({type(exc).__name__}); "
                       "used the evolutionary heuristic instead. " + result.note)
        return result


def compute_knockout_design(
    model: cobra.Model,
    product_id: str,
    *,
    method: str = "optknock",
    solver: str = "auto",
    approach: str = "best",
    biomass_id: Optional[str] = None,
    max_knockouts: int = 3,
    max_solutions: int = 5,
    time_limit: int = 120,
    min_growth_fraction: float = 0.1,
    exclude_exchanges: bool = True,
    exclude_transport: bool = False,
    exclude_blocked: bool = True,
    **_ignored,
) -> StrainDesignResult:
    """Exact MILP knockout design: OptKnock, RobustKnock or OptCouple (StrainDesign).

    ``min_growth_fraction`` forces the designed strain to keep growth at least that
    fraction of wild-type (so solutions are not lethal). Exchange/transport/blocked/
    biomass/product/ATPM reactions are protected from knockout. ``approach`` is
    'best' (optimal), 'any' (faster, suboptimal) or 'populate' (diverse set).
    """
    import straindesign as sd

    if not model.reactions.has_id(product_id):
        raise StrainDesignError(f"No product reaction '{product_id}' in the model.")
    biomass_id = biomass_id or guess_biomass_reaction(model)
    if not biomass_id or not model.reactions.has_id(biomass_id):
        raise StrainDesignError("Could not identify a valid biomass/growth reaction.")

    constraints = []
    if min_growth_fraction and min_growth_fraction > 0:
        wt = _wt_growth(model, biomass_id)
        if wt > 1e-9:
            constraints.append(f"{biomass_id} >= {min_growth_fraction * wt:.6f}")

    knockable = _knockable_reactions(
        model, product_id, biomass_id, exclude_exchanges=exclude_exchanges,
        exclude_transport=exclude_transport, exclude_blocked=exclude_blocked)
    ko_cost = {rid: 1.0 for rid in knockable}
    if not ko_cost:
        raise StrainDesignError("No reactions are available for knockout after exclusions.")

    module_type = {
        "robustknock": sd.ROBUSTKNOCK, "optcouple": sd.OPTCOUPLE,
    }.get(method, sd.OPTKNOCK)
    approach_const = {"any": sd.ANY, "populate": sd.POPULATE}.get(approach, sd.BEST)

    # Resolve solver: prefer a genuinely capable MILP solver (SCIP is bundled) over
    # GLPK, which is slow/unstable on genome-scale bi-level problems. GLPK is only a
    # last resort. (RobustKnock can hit an internal StrainDesign issue on SCIP; that
    # is caught by run_knockout_design, which falls back to the heuristic rather than
    # dropping to GLPK, which would hang.)
    chosen = solver if solver and solver != "auto" else (_capable_milp_solver() or "glpk")

    try:
        if method == "optcouple":
            # OptCouple couples a product exchange to growth: inner_objective=biomass,
            # prod_id=product (it has no separate outer objective).
            module_kwargs = dict(inner_objective=biomass_id, prod_id=product_id,
                                 constraints=constraints or None)
        else:
            module_kwargs = dict(inner_objective=biomass_id, outer_objective=product_id,
                                 constraints=constraints or None)
        with _quiet_solver():
            module = sd.SDModule(model, module_type, **module_kwargs)
            sols = sd.compute_strain_designs(
                model, sd_modules=[module], max_cost=max_knockouts,
                max_solutions=max_solutions, solution_approach=approach_const,
                solver=chosen, ko_cost=ko_cost, time_limit=time_limit)
    except Exception as exc:  # noqa: BLE001
        raise StrainDesignError(f"Strain design failed:\n{exc}") from exc

    designs = sols.get_reaction_sd() if sols is not None else []
    rows = []
    for i, design in enumerate(designs, start=1):
        knockouts = sorted(rid for rid, v in design.items() if v < 0 or v == 0)
        rows.append(_design_row(i, knockouts, model, biomass_id, product_id))
    label = method.upper()
    if not rows:
        return StrainDesignResult(
            method=label,
            note=(f"No strategy with ≤{max_knockouts} knockouts found (solver: {chosen}). "
                  "Try more knockouts, a lower minimum growth, or the heuristic method."))
    note = f"Solver: {chosen} · approach: {approach}. " + _coupling_summary(rows)
    return StrainDesignResult(method=label, table=pd.DataFrame(rows), note=note)


def heuristic_knockout_design(
    model: cobra.Model,
    product_id: str,
    *,
    biomass_id: Optional[str] = None,
    max_knockouts: int = 3,
    max_solutions: int = 5,
    min_growth_fraction: float = 0.1,
    exclude_exchanges: bool = True,
    exclude_transport: bool = False,
    exclude_blocked: bool = True,
    population: int = 40,
    generations: int = 25,
    time_limit: int = 120,
    seed: int = 0,
    **_ignored,
) -> StrainDesignResult:
    """Evolutionary (metaheuristic) knockout search.

    A fast, approximate alternative for when the exact MILP methods are intractable
    or deadlock (the role played by FastKnock / OptGene / PSO). It evolves sets of up
    to ``max_knockouts`` reactions, scoring each by the product flux attainable at the
    growth-optimal state (requiring growth ≥ a fraction of wild-type).
    """
    import random
    import time

    if not model.reactions.has_id(product_id):
        raise StrainDesignError(f"No product reaction '{product_id}' in the model.")
    biomass_id = biomass_id or guess_biomass_reaction(model)
    if not biomass_id or not model.reactions.has_id(biomass_id):
        raise StrainDesignError("Could not identify a valid biomass/growth reaction.")

    wt = _wt_growth(model, biomass_id)
    if wt <= 1e-9:
        raise StrainDesignError("Wild-type model does not grow under the current medium.")
    min_growth = min_growth_fraction * wt
    candidates = sorted(_knockable_reactions(
        model, product_id, biomass_id, exclude_exchanges=exclude_exchanges,
        exclude_transport=exclude_transport, exclude_blocked=exclude_blocked))
    if not candidates:
        raise StrainDesignError("No reactions are available for knockout after exclusions.")

    rng = random.Random(seed)
    cache: dict = {}
    work = model.copy()
    bm = work.reactions.get_by_id(biomass_id)
    prod = work.reactions.get_by_id(product_id)

    def score(kos):
        key = frozenset(kos)
        if key in cache:
            return cache[key]
        with work:
            for rid in kos:
                work.reactions.get_by_id(rid).knock_out()
            work.objective = bm
            work.objective_direction = "max"
            g = work.slim_optimize()
            if g is None or g != g or g < min_growth:
                cache[key] = (-1e9, 0.0, 0.0)
                return cache[key]
            bm.lower_bound = max(bm.lower_bound, g * 0.999)
            work.objective = prod
            work.objective_direction = "max"
            p = work.slim_optimize()
            p = float(p) if (p is not None and p == p) else 0.0
            cache[key] = (p, float(g), p)
        return cache[key]

    def random_individual():
        k = rng.randint(1, max_knockouts)
        return frozenset(rng.sample(candidates, min(k, len(candidates))))

    def mutate(ind):
        s = set(ind)
        op = rng.random()
        if op < 0.4 and len(s) < max_knockouts:
            s.add(rng.choice(candidates))
        elif op < 0.7 and s:
            s.discard(rng.choice(tuple(s)))
        elif s:
            s.discard(rng.choice(tuple(s)))
            s.add(rng.choice(candidates))
        return frozenset(s) if s else random_individual()

    def crossover(a, b):
        pool = list(a | b)
        rng.shuffle(pool)
        return frozenset(pool[:max_knockouts])

    pop = {random_individual() for _ in range(population)}
    start = time.time()
    for _ in range(generations):
        if time.time() - start > time_limit:
            break
        ranked = sorted(pop, key=lambda ind: score(ind)[0], reverse=True)
        survivors = ranked[: max(2, population // 4)]
        children = set(survivors)
        while len(children) < population and time.time() - start <= time_limit:
            a, b = rng.choice(survivors), rng.choice(survivors)
            child = mutate(crossover(a, b))
            children.add(child)
        pop = children

    scored = []
    seen = set()
    for ind in sorted(pop, key=lambda i: score(i)[0], reverse=True):
        p, g, _ = score(ind)
        if p <= 1e-6 or g < min_growth or ind in seen:
            continue
        seen.add(ind)
        scored.append((ind, g, p))
        if len(scored) >= max_solutions:
            break

    if not scored:
        return StrainDesignResult(
            method="Heuristic",
            note="No growth-coupled knockout set found. Try more knockouts, a lower "
                 "minimum growth, or a larger population/more generations.")
    rows = [_design_row(i + 1, sorted(ind), model, biomass_id, product_id)
            for i, (ind, _g, _p) in enumerate(scored)]
    note = f"Evolutionary search · {len(cache)} candidates evaluated. " + _coupling_summary(rows)
    return StrainDesignResult(method="Heuristic", table=pd.DataFrame(rows), note=note)


_COUPLING_TOL = 1e-4


def _design_row(i, knockouts, model, biomass_id, product_id) -> dict:
    """One results-table row for a knockout design, with growth-coupling honesty."""
    growth, product, coupled_min = _evaluate_design(model, knockouts, biomass_id, product_id)
    coupled = bool(coupled_min == coupled_min and coupled_min > _COUPLING_TOL)
    return {
        "solution": i,
        "knockouts": ", ".join(knockouts) if knockouts else "(none)",
        "n_knockouts": len(knockouts),
        "predicted_growth": round(growth, 6) if growth == growth else growth,
        "product_at_max_growth": round(product, 6) if product == product else product,
        "guaranteed_product": round(coupled_min, 6) if coupled_min == coupled_min else coupled_min,
        "growth_coupled": "yes" if coupled else "no",
    }


def _coupling_summary(rows: list) -> str:
    """Plain-language verdict on whether any design is meaningfully growth-coupled."""
    n_coupled = sum(1 for r in rows if r.get("growth_coupled") == "yes")
    if n_coupled:
        return (f"{n_coupled}/{len(rows)} design(s) are growth-coupled "
                f"(guaranteed product > {_COUPLING_TOL:g} at max growth).")
    return ("⚠ No design is growth-coupled of practical value (guaranteed product at max "
            f"growth ≤ {_COUPLING_TOL:g}). The product is not growth-coupled under this "
            "medium — consider two-stage/dynamic control or a different target.")


def _evaluate_design(model, knockouts, biomass_id, product_id):
    """Apply knockouts on a copy and characterise the design.

    Returns ``(growth, product_at_max_growth, coupled_min_product)`` where the
    last value is the *guaranteed* product flux at maximum growth (Issue 9): the
    minimum product flux the cell must carry while growing optimally. A design is
    only truly growth-coupled when this lower bound is meaningfully positive — a
    high ``product_at_max_growth`` alone can be an optimistic alternate optimum.
    """
    with model:
        for rid in knockouts:
            if model.reactions.has_id(rid):
                model.reactions.get_by_id(rid).knock_out()
        model.objective = biomass_id
        sol = model.optimize()
        if sol.status != "optimal":
            return float("nan"), float("nan"), float("nan")
        growth = float(sol.objective_value)
        product = float(sol.fluxes.get(product_id, float("nan")))
        # Guaranteed product at (near-)max growth: fix growth, minimise product.
        coupled_min = float("nan")
        try:
            with model:
                bm = model.reactions.get_by_id(biomass_id)
                bm.lower_bound = max(bm.lower_bound, growth - 1e-6)
                model.objective = product_id
                model.objective_direction = "min"
                mn = model.slim_optimize()
                coupled_min = float(mn) if mn is not None and mn == mn else float("nan")
        except Exception:  # noqa: BLE001
            coupled_min = float("nan")
        return growth, product, coupled_min


# --- FSEOF (over/under-expression targets) ----------------------------------
def run_fseof(
    model: cobra.Model,
    product_id: str,
    *,
    biomass_id: Optional[str] = None,
    n_steps: int = 10,
    tolerance: float = 1e-6,
    include_trending: bool = True,
    trend_correlation: float = 0.8,
    pathway_reactions: Optional[List[str]] = None,
) -> StrainDesignResult:
    """Flux Scanning with Enforced Objective Flux.

    Enforces increasing fractions of the maximum product flux while maximizing
    growth, and classifies each reaction by how its flux responds:

    * **overexpression** — flux rises monotonically as product is enforced,
    * **knockdown/deletion** — flux falls toward zero,

    which nominates amplification and down-regulation targets.

    ``tolerance`` (Issue 11) is the per-step slack allowed when testing strict
    monotonicity — exposed so near-flat real trends are not lost to a hard-coded
    value. With ``include_trending`` (default), reactions whose magnitude is not
    strictly monotone but is well correlated with the enforced flux (``|Pearson|``
    ≥ ``trend_correlation``) are still reported, marked ``… (trend)``, rather than
    silently discarded.
    """
    if not model.reactions.has_id(product_id):
        raise StrainDesignError(f"No product reaction '{product_id}' in the model.")
    biomass_id = biomass_id or guess_biomass_reaction(model)

    work = model.copy()
    # Knowing which reactions belong to the designed route lets every hit be classified
    # as pathway / precursor supply / competing sink instead of one flat list (VI.11).
    _fseof_ctx = None
    if pathway_reactions:
        _cons, _prod = _pathway_context(work, pathway_reactions)
        _fseof_ctx = (set(pathway_reactions), _cons, _prod)
    # Maximum theoretical product flux.
    work.objective = product_id
    max_prod = work.slim_optimize()
    if max_prod is None or not np.isfinite(max_prod) or abs(max_prod) < 1e-9:
        raise StrainDesignError(
            f"The model cannot produce {product_id} under the current medium.")

    fractions = np.linspace(0.1, 0.9, n_steps)
    flux_profiles = {}
    work.objective = biomass_id
    product_rxn = work.reactions.get_by_id(product_id)
    for frac in fractions:
        with work:
            enforced = frac * max_prod
            # Enforce at least this much product (respecting sign).
            if max_prod > 0:
                product_rxn.lower_bound = enforced
            else:
                product_rxn.upper_bound = enforced
            sol = work.optimize()
            if sol.status == "optimal":
                flux_profiles[frac] = sol.fluxes

    if len(flux_profiles) < 2:
        raise StrainDesignError("FSEOF could not find enough feasible enforced-flux states.")

    fdf = pd.DataFrame(flux_profiles)  # index: reactions, columns: fractions
    enforced_levels = np.array(sorted(flux_profiles.keys()), dtype=float)
    rows = []
    n_trending = 0
    n_biomass = 0
    for rid in fdf.index:
        series = fdf.loc[rid].values
        if np.all(np.abs(series) < 1e-7):
            continue
        if work.reactions.has_id(rid) and physiology.is_biomass_reaction(
                work.reactions.get_by_id(rid)):
            # Biomass draws on every precursor, so it always trends with an enforced
            # product flux — and "down-regulate biomass" is not an engineering target.
            n_biomass += 1
            continue
        start, end = series[0], series[-1]
        trend = end - start
        mags = np.abs(series)
        # Strictly monotonic (within tolerance) increase/decrease in magnitude.
        increasing = np.all(np.diff(mags) >= -tolerance) and abs(end) > abs(start) + tolerance
        decreasing = np.all(np.diff(mags) <= tolerance) and abs(end) < abs(start) - tolerance
        # Correlation of magnitude with the enforced product level (trend strength).
        corr = 0.0
        if np.std(mags) > 1e-12:
            corr = float(np.corrcoef(enforced_levels, mags)[0, 1])
        monotonic = bool(increasing or decreasing)
        if increasing:
            target = "overexpression"
        elif decreasing:
            target = "knockdown/deletion"
        elif include_trending and abs(corr) >= trend_correlation and abs(trend) > tolerance:
            # Not strictly monotone, but a strong trend — keep it, clearly marked.
            target = ("overexpression (trend)" if corr > 0 else "knockdown/deletion (trend)")
            n_trending += 1
        else:
            continue
        rows.append({
            "reaction": rid,
            "target_type": target,
            "role": _fseof_role(work, rid, _fseof_ctx),
            "monotonic": monotonic,
            "correlation": round(corr, 3),
            "flux_start": round(float(start), 4),
            "flux_end": round(float(end), 4),
            "abs_change": round(float(abs(trend)), 4),
        })
    cols = ["reaction", "target_type", "role", "monotonic", "correlation",
            "flux_start", "flux_end", "abs_change"]
    table = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
    if len(table):
        # Rank by the SIZE of the flux change, not by correlation (VI.11). Correlation
        # saturates at 1.0 for anything that merely tracks turnover — water exchange and
        # transport top the old ranking — whereas |Δflux| finds reactions that actually
        # move meaningful carbon. Roles are ordered so pathway-relevant groups come first.
        order = {"precursor supply": 0, "competing sink": 1, "pathway": 2,
                 "transport": 3, "exchange": 4, "other": 5}
        table["_r"] = table["role"].map(lambda r: order.get(r, 5))
        table = (table.sort_values(["_r", "abs_change"], ascending=[True, False])
                      .drop(columns="_r").reset_index(drop=True))
    n_generic = int((table["role"].isin(("transport", "exchange"))).sum()) if len(table) else 0
    note = (f"Scanned {len(flux_profiles)} enforced-flux levels up to {max_prod:.4g} "
            f"({product_id}). Ranked by |Δflux| within role groups.")
    if n_generic:
        note += (f" {n_generic} transport/exchange reaction(s) are listed last — they "
                 "usually track overall turnover rather than being useful targets.")
    if n_biomass:
        note += (f" {n_biomass} biomass/growth reaction(s) were excluded — they consume "
                 "every precursor, so they trend with anything and cannot be engineered.")
    if n_trending:
        note += (f" {n_trending} reaction(s) are trending but not strictly monotone "
                 f"(|corr|≥{trend_correlation}), marked '(trend)'.")
    return StrainDesignResult(method="FSEOF", table=table, note=note)


# Metabolites whose transport/exchange dominates any flux scan without being a target.
_GENERIC_MET_BASES = {"h2o", "h", "co2", "o2", "pi", "ppi", "nh4", "hco3", "photon",
                      "h2o2", "so4", "k", "na1", "cl", "mg2", "fe2", "fe3", "ca2"}


def _pathway_context(model, pathway_reactions):
    """Metabolites the designed route consumes and produces, for role classification."""
    consumed, produced = set(), set()
    for rid in (pathway_reactions or []):
        if not model.reactions.has_id(rid):
            continue
        for m, c in model.reactions.get_by_id(rid).metabolites.items():
            base = (m.id or "").lower()
            base = base.rsplit("_", 1)[0] if "_" in base else base
            if base in _GENERIC_MET_BASES:
                continue
            (consumed if c < 0 else produced).add(m.id)
    return consumed, produced


def _fseof_role(model, rid: str, context=None) -> str:
    """Classify an FSEOF hit so the output can be grouped rather than a flat list (VI.11).

    With the designed route's metabolites as ``context``:

    * **pathway** — a step of the route itself;
    * **precursor supply** — makes something the route consumes ⇒ amplify it;
    * **competing sink** — also consumes something the route needs ⇒ knock it down;
    * **transport / exchange** — water, protons, ions: statistical passengers, listed last.
    """
    try:
        rxn = model.reactions.get_by_id(rid)
    except Exception:  # noqa: BLE001
        return "other"
    if context is not None:
        route_ids, consumed, produced = context
        if rid in route_ids:
            return "pathway"
        makes_needed = consumes_needed = False
        for m, c in rxn.metabolites.items():
            if m.id in consumed:
                if c > 0:
                    makes_needed = True
                elif c < 0:
                    consumes_needed = True
        if makes_needed:
            return "precursor supply"
        if consumes_needed:
            return "competing sink"
    bases = set()
    for m in rxn.metabolites:
        mid = (m.id or "").lower()
        bases.add(mid.rsplit("_", 1)[0] if "_" in mid else mid)
    if rxn.boundary:
        return "exchange" if bases & _GENERIC_MET_BASES else "other"
    # Transport: the same species appears in more than one compartment.
    comps = {}
    for m in rxn.metabolites:
        b = (m.id or "").lower().rsplit("_", 1)[0]
        comps.setdefault(b, set()).add(getattr(m, "compartment", "") or "")
    if any(len(c) > 1 for c in comps.values()):
        return "transport" if bases & _GENERIC_MET_BASES else "other"
    if bases <= _GENERIC_MET_BASES:
        return "other"
    return "pathway"


def group_fseof_targets(result, *, hide_generic: bool = True):
    """Split an FSEOF table into the groups the user should read, in order (VI.11).

    Returns ``{"precursor supply": df, "competing sink": df, …}``. With
    ``hide_generic`` the transport/exchange passengers are omitted — the GUI exposes this
    as a "show all" toggle so nothing is hidden irrecoverably.
    """
    import pandas as pd
    table = getattr(result, "table", None)
    if table is None or not len(table) or "role" not in table.columns:
        return {}
    groups = {}
    for role in ("pathway", "precursor supply", "competing sink", "transport",
                 "exchange", "other"):
        sub = table[table["role"] == role]
        if not len(sub):
            continue
        if hide_generic and role in ("transport", "exchange"):
            continue
        groups[role] = sub.reset_index(drop=True)
    return groups


def run_metabolite_overproduction(
    model: cobra.Model,
    metabolite_id: str,
    *,
    biomass_id: Optional[str] = None,
    n_steps: int = 10,
    tolerance: float = 1e-6,
    include_trending: bool = True,
    secretion: bool = False,
    gas: bool = False,
) -> StrainDesignResult:
    """Find the reactions (and directions) to engineer for accumulating a *metabolite*.

    Many designs target a product *reaction*, but the real goal is usually to
    accumulate a *metabolite*. By default this adds a demand reaction that drains
    the target (forcing its net production). With ``secretion`` (Issue 10) it
    instead builds a realistic export route — transport to the extracellular
    compartment plus an exchange (a gas-phase sink when ``gas``) — so the enforced
    product reflects real secretion rather than an internal drain. FSEOF is then
    run on that reaction to nominate amplification / knock-down targets.
    """
    if not model.metabolites.has_id(metabolite_id):
        raise StrainDesignError(f"No metabolite '{metabolite_id}' in the model.")

    work = model.copy()
    met = work.metabolites.get_by_id(metabolite_id)
    if secretion:
        product_id = secrete_metabolite(work, met.id, gas=gas)
    else:
        demand_id = f"DM_{met.id}"
        if work.reactions.has_id(demand_id):
            work.reactions.get_by_id(demand_id).bounds = (0.0, 1000.0)
            product_id = demand_id
        else:
            product_id = work.add_boundary(met, type="demand").id  # met --> (net production)

    result = run_fseof(work, product_id, biomass_id=biomass_id, n_steps=n_steps,
                       tolerance=tolerance, include_trending=include_trending)

    # Annotate each nominated reaction with the direction it should run.
    if not result.table.empty:
        result.table.insert(
            2, "direction",
            result.table["flux_end"].apply(
                lambda f: "forward" if f > 1e-9 else ("reverse" if f < -1e-9 else "—")))
        result.table = result.table[result.table["reaction"] != product_id].reset_index(drop=True)
    result.method = "Metabolite overproduction (FSEOF)"
    route = ("secretion (transport + exchange"
             + (" + gas sink" if gas else "") + ")") if secretion else "demand drain"
    result.note = (f"Targets to accumulate {metabolite_id} via {route}. " + result.note +
                   " 'direction' shows which way each reaction must run; 'overexpression' "
                   "reactions should be amplified, 'knockdown/deletion' reduced.")
    return result
