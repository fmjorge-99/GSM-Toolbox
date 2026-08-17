"""Multi-organism (community) metabolic modeling (Phase 6).

Combines several single-organism models into one community model using the
standard compartmentalized approach:

* each organism's *intracellular* reactions and metabolites are tagged with the
  organism name (so they stay private),
* *extracellular* metabolites (those involved in exchange reactions) are pooled
  into a single shared environment that all organisms draw from and secrete into,
* one community-level exchange reaction per shared metabolite connects the pool to
  the environment,
* the objective becomes the sum of the members' biomass reactions.

Community FBA then predicts how the organisms grow together and trade metabolites.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import cobra
import pandas as pd

from .editing import guess_biomass_reaction


class CommunityError(Exception):
    """Raised when a community model cannot be built."""


@dataclass
class CommunityModel:
    model: cobra.Model
    member_names: List[str]
    member_biomass: Dict[str, str] = field(default_factory=dict)  # name -> tagged biomass id


def _external_metabolite_ids(model: cobra.Model) -> set:
    """External metabolites = those consumed/produced by boundary (exchange) reactions."""
    ext = set()
    for rxn in model.reactions:
        if rxn.boundary:
            ext.update(m.id for m in rxn.metabolites)
    return ext


def build_community(models: List[cobra.Model], names: Optional[List[str]] = None) -> CommunityModel:
    """Combine ``models`` into one community model with a shared environment."""
    if len(models) < 2:
        raise CommunityError("A community needs at least two models.")
    if names is None:
        names = [(m.id or f"org{i+1}") for i, m in enumerate(models)]
    # Ensure unique, file-system-friendly tags.
    seen = {}
    tags = []
    for n in names:
        base = (n or "org").replace(" ", "_")
        k = seen.get(base, 0)
        seen[base] = k + 1
        tags.append(base if k == 0 else f"{base}{k+1}")

    comm = cobra.Model("community")
    shared: Dict[str, cobra.Metabolite] = {}
    ex_bounds: Dict[str, Tuple[float, float]] = {}
    member_biomass: Dict[str, str] = {}

    for model, tag in zip(models, tags):
        ext_ids = _external_metabolite_ids(model)
        biomass = guess_biomass_reaction(model)
        local_map: Dict[str, cobra.Metabolite] = {}

        for met in model.metabolites:
            if met.id in ext_ids:
                if met.id not in shared:
                    shared[met.id] = cobra.Metabolite(
                        met.id, name=met.name, formula=met.formula,
                        charge=met.charge, compartment=met.compartment or "e")
                local_map[met.id] = shared[met.id]
            else:
                tid = f"{met.id}__{tag}"
                local_map[met.id] = cobra.Metabolite(
                    tid, name=f"{met.name} [{tag}]", formula=met.formula,
                    charge=met.charge, compartment=f"{met.compartment or 'c'}__{tag}")

        new_reactions = []
        for rxn in model.reactions:
            if rxn.boundary:
                # Accumulate the most permissive exchange bounds for the shared pool.
                for met in rxn.metabolites:
                    lb, ub = ex_bounds.get(met.id, (0.0, 0.0))
                    ex_bounds[met.id] = (min(lb, rxn.lower_bound), max(ub, rxn.upper_bound))
                continue
            nr = cobra.Reaction(f"{rxn.id}__{tag}")
            nr.name = f"{rxn.name} [{tag}]"
            nr.bounds = (rxn.lower_bound, rxn.upper_bound)
            new_reactions.append((nr, rxn))

        comm.add_reactions([nr for nr, _ in new_reactions])
        for nr, rxn in new_reactions:
            nr.add_metabolites({local_map[m.id]: c for m, c in rxn.metabolites.items()})

        if biomass:
            tagged_bm = f"{biomass}__{tag}"
            if comm.reactions.has_id(tagged_bm):
                member_biomass[tag] = tagged_bm

    # Community exchanges on the shared pool.
    for mid, met in shared.items():
        ex = cobra.Reaction(f"EX_{mid}")
        lb, ub = ex_bounds.get(mid, (0.0, 1000.0))
        ex.bounds = (lb, ub if ub != 0 else 1000.0)
        comm.add_reactions([ex])
        ex.add_metabolites({met: -1})

    if not member_biomass:
        raise CommunityError("Could not identify biomass reactions in the member models.")
    comm.objective = {comm.reactions.get_by_id(b): 1.0 for b in member_biomass.values()}
    comm.id = "community_" + "_".join(tags)
    return CommunityModel(model=comm, member_names=tags, member_biomass=member_biomass)


def member_solo_max(model: cobra.Model, member_biomass: Dict[str, str]) -> Dict[str, float]:
    """The maximum growth each member could reach *in the community context* (i.e.
    if it alone were optimised while the others are free). Used as the reference for
    the per-member minimum-growth floor, so 'grow at least X% of what you could'."""
    solo: Dict[str, float] = {}
    for name, bm in member_biomass.items():
        if not model.reactions.has_id(bm):
            solo[name] = 0.0
            continue
        with model:
            model.objective = bm
            s = model.optimize()
            solo[name] = float(s.objective_value or 0.0) if s.status == "optimal" else 0.0
    return solo


def _apply_community_objective(model, member_biomass, weights, min_growth_fraction):
    """Set a *weighted* biomass objective (dominance) and optional per-member growth
    floors (no member is starved). Returns (normalised_weights, floors)."""
    names = list(member_biomass)
    w = {n: 1.0 for n in names}
    if weights:
        for k, v in weights.items():
            if k in w:
                w[k] = max(0.0, float(v))
    if not any(w.values()):                      # guard against all-zero weights
        w = {n: 1.0 for n in names}
    total_w = sum(w.values()) or 1.0
    w = {k: v / total_w * len(names) for k, v in w.items()}   # normalise to mean 1.0

    floors: Dict[str, float] = {}
    if min_growth_fraction and min_growth_fraction > 0:
        solo = member_solo_max(model, member_biomass)
        for n, bm in member_biomass.items():
            if model.reactions.has_id(bm):
                r = model.reactions.get_by_id(bm)
                floor = float(min_growth_fraction) * solo.get(n, 0.0)
                if floor > r.lower_bound:
                    r.lower_bound = floor
                floors[n] = floor
    model.objective = {model.reactions.get_by_id(bm): w[n]
                       for n, bm in member_biomass.items() if model.reactions.has_id(bm)}
    return w, floors


def member_growth_table(model: cobra.Model, member_biomass: Dict[str, str], *,
                        weights: Optional[Dict[str, float]] = None,
                        min_growth_fraction: float = 0.0) -> pd.DataFrame:
    """Optimize ``model`` and report each member's biomass flux + the total.

    ``weights`` set each member's relative dominance in the community objective;
    ``min_growth_fraction`` forces every member to grow at >= that fraction of its
    own community-context maximum (so the faster grower can't starve the others)."""
    with model:
        w, floors = _apply_community_objective(model, member_biomass, weights, min_growth_fraction)
        sol = model.optimize()
        feasible = sol.status == "optimal"
        rows = []
        total = 0.0
        for name, bm_id in member_biomass.items():
            flux = float(sol.fluxes.get(bm_id, 0.0)) if feasible else float("nan")
            total += flux if flux == flux else 0.0
            rows.append({"member": name, "biomass_reaction": bm_id, "growth": flux,
                         "weight": round(w.get(name, 1.0), 3),
                         "min_growth": round(floors.get(name, 0.0), 5)})
        rows.append({"member": "TOTAL (community)", "biomass_reaction": "", "growth": total,
                     "weight": "", "min_growth": ""})
    df = pd.DataFrame(rows)
    if not feasible:
        df.attrs["infeasible"] = (
            "No feasible community state at these settings — the minimum-growth "
            "requirement is too high for the shared nutrients. Lower it and re-run.")
    return df


def community_fba(community: CommunityModel, *, weights: Optional[Dict[str, float]] = None,
                  min_growth_fraction: float = 0.0) -> pd.DataFrame:
    """Run FBA on the community and report each member's growth + the total."""
    return member_growth_table(community.model, community.member_biomass,
                               weights=weights, min_growth_fraction=min_growth_fraction)
