"""Reaction / metabolite / objective editing operations on a cobra.Model.

These functions are deliberately *pure operations on a model*: they mutate the
model they are given and raise :class:`EditError` on invalid input. Undo/redo
and modified-item tracking live in :class:`gsm_toolbox.core.project.Project`,
which snapshots the model around each edit.
"""

from __future__ import annotations

from typing import Dict, Optional

import cobra


class EditError(Exception):
    """Raised when a requested edit is invalid."""


def add_reaction(
    model: cobra.Model,
    reaction_id: str,
    reaction_string: str,
    *,
    name: str = "",
    lower_bound: float = 0.0,
    upper_bound: float = 1000.0,
    gene_reaction_rule: str = "",
    subsystem: str = "",
) -> cobra.Reaction:
    """Add a new reaction defined by a human-readable reaction string.

    ``reaction_string`` uses COBRApy's :meth:`Reaction.build_reaction_from_string`
    syntax, e.g. ``"atp_c + h2o_c --> adp_c + pi_c + h_c"``. Metabolites that do
    not yet exist are created automatically (assigned to the first compartment).
    """
    reaction_id = reaction_id.strip()
    if not reaction_id:
        raise EditError("Reaction ID cannot be empty.")
    if model.reactions.has_id(reaction_id):
        raise EditError(f"A reaction with ID '{reaction_id}' already exists.")
    if lower_bound > upper_bound:
        raise EditError("Lower bound cannot be greater than upper bound.")

    rxn = cobra.Reaction(reaction_id)
    rxn.name = name
    rxn.subsystem = subsystem
    rxn.lower_bound = lower_bound
    rxn.upper_bound = upper_bound
    model.add_reactions([rxn])

    try:
        if reaction_string.strip():
            rxn.build_reaction_from_string(reaction_string)
    except Exception as exc:  # noqa: BLE001
        model.remove_reactions([rxn])
        raise EditError(f"Could not parse the reaction equation:\n{exc}") from exc

    if gene_reaction_rule.strip():
        try:
            rxn.gene_reaction_rule = gene_reaction_rule.strip()
        except Exception as exc:  # noqa: BLE001
            raise EditError(f"Invalid gene-reaction rule:\n{exc}") from exc
    return rxn


def remove_reaction(model: cobra.Model, reaction_id: str, *, remove_orphans: bool = True) -> None:
    """Remove a reaction by ID; optionally drop now-orphaned metabolites/genes."""
    if not model.reactions.has_id(reaction_id):
        raise EditError(f"No reaction with ID '{reaction_id}'.")
    rxn = model.reactions.get_by_id(reaction_id)
    model.remove_reactions([rxn], remove_orphans=remove_orphans)


def set_bounds(model: cobra.Model, reaction_id: str, lower: float, upper: float) -> None:
    """Set the flux bounds of a reaction."""
    if lower > upper:
        raise EditError("Lower bound cannot be greater than upper bound.")
    if not model.reactions.has_id(reaction_id):
        raise EditError(f"No reaction with ID '{reaction_id}'.")
    model.reactions.get_by_id(reaction_id).bounds = (lower, upper)


def set_gene_reaction_rule(model: cobra.Model, reaction_id: str, rule: str) -> None:
    """Update the gene-reaction rule (GPR) of a reaction."""
    if not model.reactions.has_id(reaction_id):
        raise EditError(f"No reaction with ID '{reaction_id}'.")
    try:
        model.reactions.get_by_id(reaction_id).gene_reaction_rule = rule.strip()
    except Exception as exc:  # noqa: BLE001
        raise EditError(f"Invalid gene-reaction rule:\n{exc}") from exc


def set_objective(model: cobra.Model, reaction_id: str, direction: str = "max") -> None:
    """Set the model objective to a single reaction and the optimization sense."""
    if not model.reactions.has_id(reaction_id):
        raise EditError(f"No reaction with ID '{reaction_id}'.")
    if direction not in ("max", "min"):
        raise EditError("Direction must be 'max' or 'min'.")
    model.objective = reaction_id
    model.objective_direction = direction


def set_weighted_objective(
    model: cobra.Model, terms: Dict[str, float], direction: str = "max"
) -> None:
    """Set a (possibly multi-reaction) linear objective.

    ``terms`` maps reaction id -> weight, e.g. ``{"BIOMASS": 1.0, "EX_succ_e": 0.5}``
    to balance growth with succinate production. Weights of 0 are ignored.
    """
    if direction not in ("max", "min"):
        raise EditError("Direction must be 'max' or 'min'.")
    active = {rid: w for rid, w in terms.items() if w}
    if not active:
        raise EditError("Provide at least one objective reaction with a non-zero weight.")
    obj = {}
    for rid, weight in active.items():
        if not model.reactions.has_id(rid):
            raise EditError(f"No reaction with ID '{rid}'.")
        obj[model.reactions.get_by_id(rid)] = float(weight)
    model.objective = obj
    model.objective_direction = direction


def current_objective_terms(model: cobra.Model) -> Dict[str, float]:
    """Return the model's current objective as ``{reaction_id: coefficient}``."""
    return {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient}


def guess_biomass_reaction(model: cobra.Model) -> Optional[str]:
    """Best-effort detection of the biomass/growth reaction.

    Prefers a name/id match for "biomass" so the result is stable even after the
    user changes the model objective to a product. Falls back to the current
    objective reaction when no biomass-named reaction exists.
    """
    for rxn in model.reactions:
        if "biomass" in rxn.id.lower() or "biomass" in (rxn.name or "").lower():
            return rxn.id
    current = current_objective_terms(model)
    if current:
        return max(current, key=lambda k: abs(current[k]))
    return None


def max_abs_flux(model: cobra.Model, reaction_id: str) -> float:
    """Largest attainable |flux| of a reaction (maximize then minimize it).

    Used to normalize weighted objectives so percentage weights are meaningful:
    biomass (~1) and a product exchange (~20) can then be combined on a comparable
    scale. Returns 0.0 if the reaction cannot carry flux.
    """
    if not model.reactions.has_id(reaction_id):
        return 0.0
    best = 0.0
    rxn = model.reactions.get_by_id(reaction_id)
    for sense in ("max", "min"):
        with model:
            model.objective = rxn
            model.objective_direction = sense
            val = model.slim_optimize()
        if val is not None and val == val:  # not NaN
            best = max(best, abs(float(val)))
    return best


def normalized_objective_coefficients(
    model: cobra.Model, weights: Dict[str, float]
) -> Dict[str, float]:
    """Convert relative weights (e.g. 0.6/0.4) into objective coefficients.

    Each weight is divided by the reaction's max attainable flux so that the terms
    contribute in proportion to the weights rather than to their raw flux scales.
    """
    coeffs = {}
    for rid, w in weights.items():
        if not w:
            continue
        m = max_abs_flux(model, rid)
        coeffs[rid] = (w / m) if m > 1e-9 else w
    return coeffs


def reaction_to_dict(rxn: cobra.Reaction, flux: Optional[float] = None) -> Dict:
    """Flatten a reaction into a dict suitable for table display."""
    return {
        "id": rxn.id,
        "name": rxn.name or "",
        "reaction": rxn.build_reaction_string(use_metabolite_names=False),
        "lower_bound": rxn.lower_bound,
        "upper_bound": rxn.upper_bound,
        "subsystem": rxn.subsystem or "",
        "gpr": rxn.gene_reaction_rule or "",
        "flux": flux,
    }
