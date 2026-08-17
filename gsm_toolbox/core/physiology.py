"""Growth-physiology helpers: what is the model actually set up to do?

Constraint-based models encode the growth condition entirely in exchange bounds
and the objective. That makes it easy to analyse the *wrong* physiology without
noticing — e.g. running a cyanobacterial model in the dark on glucose because
that is how the SBML happened to ship (Issue 5/6). This module summarises the
active medium/objective in plain language, detects phototroph-capable models,
and builds one-click growth presets (photoautotrophic / heterotrophic).

Pure-Python, no Qt — usable headlessly and from tests.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import cobra

# clean_label lives in network_graph (a lower-level module) so display cleanup is
# shared by the tables/equations too; re-exported here for backward compatibility.
from .network_graph import clean_label  # noqa: F401

# Metabolite base ids treated as *inorganic* carbon (not an organic C source).
_INORGANIC_C = {"co2", "hco3", "co3", "cot"}
# Photon / light exchange base ids across conventions.
_PHOTON_BASES = {"photon", "photonVis", "hnu", "light"}


def _base(mid: str) -> str:
    return mid.rsplit("_", 1)[0] if "_" in mid else mid


def _carbon_count(met: cobra.Metabolite) -> int:
    total = 0
    for elem, num in re.findall(r"([A-Z][a-z]?)(\d*)", getattr(met, "formula", "") or ""):
        if elem == "C":
            total += int(num) if num else 1
    return total


def photon_exchange(model: cobra.Model) -> Optional[str]:
    """Best-effort id of a light/photon exchange reaction, if the model has one."""
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            if _base(met.id).lower() in {b.lower() for b in _PHOTON_BASES}:
                return rxn.id
        if "photon" in rxn.id.lower() or rxn.id.lower().endswith("_hnu_e"):
            return rxn.id
    return None


def active_uptakes(model: cobra.Model) -> List[dict]:
    """Exchange reactions currently allowing uptake (lower_bound < 0)."""
    rows = []
    for rxn in model.exchanges:
        if rxn.lower_bound < 0:
            met = next(iter(rxn.metabolites), None)
            rows.append({
                "id": rxn.id,
                "name": clean_label(rxn.name or (met.name if met else "") or rxn.id),
                "uptake": -rxn.lower_bound,
                "carbon": _carbon_count(met) if met else 0,
                "base": _base(met.id).lower() if met else "",
            })
    return rows


def carbon_sources(model: cobra.Model) -> List[dict]:
    """Active uptake exchanges that supply *organic* carbon (excludes CO2/HCO3)."""
    return [u for u in active_uptakes(model)
            if u["carbon"] > 0 and u["base"] not in _INORGANIC_C]


def inorganic_carbon_uptakes(model: cobra.Model) -> List[dict]:
    return [u for u in active_uptakes(model) if u["base"] in _INORGANIC_C]


def objective_reactions(model: cobra.Model) -> List[str]:
    return [r.id for r in model.reactions if r.objective_coefficient]


def is_biomass_reaction(rxn) -> bool:
    """True for a biomass/growth pseudo-reaction.

    These are never engineering targets. Biomass consumes essentially every precursor in
    the cell, so it correlates with any flux scan and competes for every intermediate —
    it tops FSEOF rankings and fills the competing-reaction list while telling the user
    nothing they can act on. "Down-regulate biomass" is not a strain design.
    """
    text = f"{getattr(rxn, 'id', '')} {getattr(rxn, 'name', '') or ''}".lower()
    if "biomass" in text or "growth" in text:
        return True
    rid = (getattr(rxn, "id", "") or "").lower()
    # `bio1`, `BIOMASS_Ec_iJO1366_core_53p95M`, `Growth` … but not `biotin` or `bioA`.
    return bool(re.match(r"^bio\d", rid))


def list_biomass_reactions(model: cobra.Model) -> List[dict]:
    """Candidate biomass/growth reactions, with a guess at auto- vs heterotrophic.

    A reaction is a biomass candidate if its id/name mentions biomass/growth. The
    'kind' hint keys off common naming (auto/photo vs hetero) so a phototroph
    model that ships several biomass reactions can be told apart."""
    rows = []
    for rxn in model.reactions:
        text = f"{rxn.id} {rxn.name or ''}".lower()
        if "biomass" in text or "growth" in text or rxn.id.lower().startswith("bio"):
            kind = ""
            if any(k in text for k in ("auto", "photo")):
                kind = "autotrophic"
            elif any(k in text for k in ("hetero", "dark", "gluc")):
                kind = "heterotrophic"
            rows.append({"id": rxn.id, "name": rxn.name or "", "kind": kind,
                         "objective": bool(rxn.objective_coefficient)})
    return rows


def is_phototroph_capable(model: cobra.Model) -> bool:
    """True when the model offers a photon exchange (so photoautotrophy is possible)."""
    return photon_exchange(model) is not None


@dataclass
class PhysiologySummary:
    objective: str = ""
    energy_source: str = "unknown"     # light | organic | inorganic/unclear
    carbon_source: str = "none"
    photon_active: bool = False
    warnings: List[str] = field(default_factory=list)
    biomass_candidates: List[dict] = field(default_factory=list)


def summarize(model: cobra.Model) -> PhysiologySummary:
    """Plain-language summary of the active objective, carbon and energy source,
    with warnings when the configuration looks physiologically inconsistent."""
    obj = " + ".join(objective_reactions(model)) or "(none)"
    photon_id = photon_exchange(model)
    photon_active = bool(photon_id and model.reactions.get_by_id(photon_id).lower_bound < 0)
    organic = carbon_sources(model)
    inorganic = inorganic_carbon_uptakes(model)

    if photon_active:
        energy = "light"
    elif organic:
        energy = "organic"
    else:
        energy = "inorganic/unclear"

    if organic:
        carbon = ", ".join(sorted({u["name"] or u["id"] for u in organic}))
    elif inorganic:
        carbon = ", ".join(sorted({u["name"] or u["id"] for u in inorganic})) + " (inorganic)"
    else:
        carbon = "none active"

    warnings = []
    # A phototroph capable of light but running in the dark on organic carbon.
    if photon_id and not photon_active and is_phototroph_capable(model) and organic:
        warnings.append("This model has a photon (light) exchange but light uptake is OFF "
                        "and an organic carbon source is active — you are analysing "
                        "heterotrophic (dark) growth. Use Growth settings for a "
                        "photoautotrophic preset if that is not intended.")
    if photon_active and organic:
        warnings.append("Both light AND an organic carbon source are active "
                        "(mixotrophic) — confirm this is intended.")
    if not organic and not inorganic and not photon_active:
        warnings.append("No carbon or energy source appears active — growth may be "
                        "infeasible. Check the medium in Growth settings.")

    return PhysiologySummary(
        objective=obj, energy_source=energy, carbon_source=carbon,
        photon_active=photon_active, warnings=warnings,
        biomass_candidates=list_biomass_reactions(model))


def apply_photoautotrophic_preset(model: cobra.Model, *, photon_uptake: float = 1000.0,
                                  co2_uptake: float = 1000.0,
                                  biomass_id: Optional[str] = None) -> List[str]:
    """Open light + inorganic carbon, close organic carbon, set an autotrophic
    biomass objective. Returns human-readable notes about what changed."""
    notes = []
    photon_id = photon_exchange(model)
    if photon_id:
        model.reactions.get_by_id(photon_id).lower_bound = -abs(photon_uptake)
        notes.append(f"Opened light uptake ({photon_id}).")
    else:
        notes.append("No photon exchange found — cannot enable light.")
    # Enable inorganic carbon (CO2 / bicarbonate).
    opened_c = []
    for rxn in model.exchanges:
        met = next(iter(rxn.metabolites), None)
        if met and _base(met.id).lower() in _INORGANIC_C:
            rxn.lower_bound = -abs(co2_uptake)
            opened_c.append(rxn.id)
    if opened_c:
        notes.append("Enabled inorganic carbon uptake: " + ", ".join(sorted(opened_c)) + ".")
    # Close organic carbon uptake.
    closed = []
    for u in carbon_sources(model):
        model.reactions.get_by_id(u["id"]).lower_bound = 0.0
        closed.append(u["id"])
    if closed:
        notes.append("Closed organic carbon uptake: " + ", ".join(sorted(closed)) + ".")
    # Set an autotrophic biomass objective if identifiable.
    if biomass_id is None:
        autos = [b for b in list_biomass_reactions(model) if b["kind"] == "autotrophic"]
        biomass_id = autos[0]["id"] if autos else None
    if biomass_id and model.reactions.has_id(biomass_id):
        model.objective = biomass_id
        notes.append(f"Set objective to autotrophic biomass ({biomass_id}).")
    return notes


def has_organic_carbon_exchange(model: cobra.Model) -> bool:
    """Whether the model has any exchange for an organic (non-CO2/HCO3) carbon source."""
    for rxn in model.exchanges:
        met = next(iter(rxn.metabolites), None)
        if met and _carbon_count(met) > 0 and _base(met.id).lower() not in _INORGANIC_C:
            return True
    return False


def available_growth_modes(model: cobra.Model) -> Dict[str, bool]:
    """Which growth modes make sense for this model (Issue R4).

    * ``autotrophic``  — needs a photon exchange (light-driven fixation possible).
    * ``mixotrophic``  — needs both light AND an organic carbon source.
    * ``heterotrophic``— needs an organic carbon source.

    A model with no photon exchange (a plain heterotroph) therefore offers only
    ``heterotrophic``; a cyanobacterial model offers all three.
    """
    photo = is_phototroph_capable(model)
    organic = has_organic_carbon_exchange(model)
    return {
        "autotrophic": photo,
        "mixotrophic": photo and organic,
        "heterotrophic": organic,
    }


def substrate_exchanges(model: cobra.Model) -> List[dict]:
    """Every exchange that could serve as a nutrient/substrate, for the picker UI.

    Returns rows ``{id, name, carbon, active}`` where ``active`` means uptake is
    currently enabled (lower_bound < 0)."""
    rows = []
    for rxn in model.exchanges:
        met = next(iter(rxn.metabolites), None)
        rows.append({
            "id": rxn.id,
            "name": clean_label(rxn.name or (met.name if met else "") or rxn.id),
            "carbon": _carbon_count(met) if met else 0,
            "active": rxn.lower_bound < 0,
        })
    rows.sort(key=lambda r: (not r["active"], r["name"].lower()))
    return rows


def apply_mixotrophic_preset(model: cobra.Model, *, photon_uptake: float = 1000.0,
                             co2_uptake: float = 1000.0,
                             biomass_id: Optional[str] = None) -> List[str]:
    """Open light AND keep the current organic carbon source(s) active, set a
    mixotrophic biomass objective when one is identifiable. Returns notes."""
    notes = []
    photon_id = photon_exchange(model)
    if photon_id:
        model.reactions.get_by_id(photon_id).lower_bound = -abs(photon_uptake)
        notes.append(f"Opened light uptake ({photon_id}).")
    for rxn in model.exchanges:
        met = next(iter(rxn.metabolites), None)
        if met and _base(met.id).lower() in _INORGANIC_C:
            rxn.lower_bound = -abs(co2_uptake)
    if not carbon_sources(model):
        notes.append("No organic carbon source is currently active — enable one in the "
                     "medium (Select Substrates) for true mixotrophy.")
    else:
        notes.append("Kept the active organic carbon source(s) alongside light.")
    if biomass_id is None:
        mixo = [b for b in list_biomass_reactions(model)
                if "mixo" in (b["id"] + b["name"]).lower()]
        biomass_id = mixo[0]["id"] if mixo else None
    if biomass_id and model.reactions.has_id(biomass_id):
        model.objective = biomass_id
        notes.append(f"Set objective to mixotrophic biomass ({biomass_id}).")
    return notes


def apply_heterotrophic_mode(model: cobra.Model, *, biomass_id: Optional[str] = None) -> List[str]:
    """Switch to heterotrophic mode: close light and set a heterotrophic biomass
    objective, leaving the organic carbon source(s) as chosen in the medium/picker.
    Returns notes."""
    notes = []
    photon_id = photon_exchange(model)
    if photon_id and model.reactions.get_by_id(photon_id).lower_bound < 0:
        model.reactions.get_by_id(photon_id).lower_bound = 0.0
        notes.append("Closed light uptake.")
    if not carbon_sources(model):
        notes.append("No organic carbon source is active — enable one via the medium "
                     "(Select Substrates) so the cell can grow heterotrophically.")
    if biomass_id is None:
        hets = [b for b in list_biomass_reactions(model) if b["kind"] == "heterotrophic"]
        biomass_id = hets[0]["id"] if hets else None
    if biomass_id and model.reactions.has_id(biomass_id):
        model.objective = biomass_id
        notes.append(f"Set objective to heterotrophic biomass ({biomass_id}).")
    if not notes:
        notes.append("Heterotrophic mode (light off).")
    return notes


def find_exchange_for_base(model: cobra.Model, base: str) -> Optional[str]:
    """Find the exchange reaction for a metabolite with the given base id (e.g.
    'glc__D', 'ac', 'glyc'), preferring the extracellular one."""
    base = base.lower()
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            if _base(met.id).lower() == base:
                return rxn.id
    return None


# Common heterotrophic carbon sources: label -> candidate metabolite base ids.
HETEROTROPH_CARBON = {
    "Glucose": ["glc__D", "glc_D", "glc"],
    "Acetate": ["ac"],
    "Glycerol": ["glyc"],
    "Succinate": ["succ"],
}


def apply_heterotrophic_preset(model: cobra.Model, carbon_label: str, *,
                               uptake: float = 10.0, aerobic: bool = True,
                               biomass_id: Optional[str] = None) -> List[str]:
    """Close light, switch to a single organic carbon source and (optionally) a
    heterotrophic biomass objective. Returns human-readable notes."""
    notes = []
    photon_id = photon_exchange(model)
    if photon_id:
        model.reactions.get_by_id(photon_id).lower_bound = 0.0
        notes.append("Closed light uptake.")
    # Close all current organic carbon uptakes first.
    for u in carbon_sources(model):
        model.reactions.get_by_id(u["id"]).lower_bound = 0.0
    # Open the chosen carbon source.
    ex_id = None
    for base in HETEROTROPH_CARBON.get(carbon_label, [carbon_label.lower()]):
        ex_id = find_exchange_for_base(model, base)
        if ex_id:
            break
    if ex_id:
        model.reactions.get_by_id(ex_id).lower_bound = -abs(uptake)
        notes.append(f"Opened {carbon_label} uptake ({ex_id}) at {uptake:g} mmol/gDW/h.")
    else:
        notes.append(f"No exchange found for {carbon_label} — carbon source unchanged.")
    # Oxygen.
    from .media import set_aerobic
    if set_aerobic(model, aerobic):
        notes.append("Aerobic." if aerobic else "Anaerobic.")
    # Heterotrophic biomass objective if identifiable.
    if biomass_id is None:
        hets = [b for b in list_biomass_reactions(model) if b["kind"] == "heterotrophic"]
        biomass_id = hets[0]["id"] if hets else None
    if biomass_id and model.reactions.has_id(biomass_id):
        model.objective = biomass_id
        notes.append(f"Set objective to heterotrophic biomass ({biomass_id}).")
    return notes
