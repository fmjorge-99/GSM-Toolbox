"""Growth medium and exchange-reaction management.

In constraint-based modeling the "medium" is defined by the uptake bounds of the
exchange reactions. COBRApy exposes this via ``model.medium`` (a dict of
``exchange_id -> max uptake rate``). This module wraps that with helpers the GUI
can use, plus aerobic/anaerobic presets.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import cobra

# Common identifiers for the oxygen exchange across model conventions.
_O2_EXCHANGE_CANDIDATES = ("EX_o2_e", "EX_o2_e_", "EX_o2(e)", "EX_o2_LPAREN_e_RPAREN_")


def list_exchanges(model: cobra.Model) -> List[dict]:
    """Return a row per exchange reaction with its current uptake/secretion bounds."""
    rows = []
    medium = model.medium
    for rxn in model.exchanges:
        rows.append(
            {
                "id": rxn.id,
                "name": rxn.name or "",
                "lower_bound": rxn.lower_bound,
                "upper_bound": rxn.upper_bound,
                # max uptake from model.medium (0 if not in medium)
                "uptake": medium.get(rxn.id, 0.0),
            }
        )
    return rows


def get_medium(model: cobra.Model) -> Dict[str, float]:
    """Return the current medium as ``{exchange_id: max_uptake}``."""
    return dict(model.medium)


def set_medium(model: cobra.Model, medium: Dict[str, float]) -> None:
    """Replace the model medium. Keys must be valid exchange reaction IDs."""
    model.medium = medium


def set_exchange_bounds(
    model: cobra.Model, exchange_id: str, lower: float, upper: float
) -> None:
    """Directly set the bounds of a single exchange reaction."""
    if not model.reactions.has_id(exchange_id):
        raise KeyError(f"No exchange reaction '{exchange_id}'.")
    model.reactions.get_by_id(exchange_id).bounds = (lower, upper)


def find_oxygen_exchange(model: cobra.Model) -> Optional[str]:
    """Best-effort lookup of the oxygen exchange reaction ID."""
    for candidate in _O2_EXCHANGE_CANDIDATES:
        if model.reactions.has_id(candidate):
            return candidate
    # Fall back to a metabolite-name based search.
    for rxn in model.exchanges:
        if any(met.id.startswith("o2_") or met.id == "o2" for met in rxn.metabolites):
            return rxn.id
    return None


def set_aerobic(model: cobra.Model, aerobic: bool, uptake: float = 1000.0) -> bool:
    """Toggle oxygen availability. Returns True if an O2 exchange was found."""
    o2 = find_oxygen_exchange(model)
    if o2 is None:
        return False
    rxn = model.reactions.get_by_id(o2)
    # Uptake is the (negative) lower bound for an exchange reaction.
    rxn.lower_bound = -abs(uptake) if aerobic else 0.0
    return True
