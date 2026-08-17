"""Reaction/metabolite categories and subset (pathway) extraction.

A :class:`Category` is a named, colored group of reactions the user defines to
explore a pathway or module of interest. Categories let the user:

* visualize a subset of the network in isolation,
* run analyses restricted to that subset.

For subset analysis we build a standalone sub-model containing only the category
reactions. Every metabolite that the category shares with the *rest* of the
network (a "cut" metabolite) is given a free exchange reaction, so the subset can
be simulated on its own with those connections treated as default input/output
fluxes — exactly the behavior requested for pathway-level analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set

import cobra

# A palette of distinct, pleasant colors assigned to new categories in order.
DEFAULT_PALETTE = [
    "#E8453C", "#1A73E8", "#34A853", "#FBBC04", "#9C27B0",
    "#00ACC1", "#FF7043", "#5E35B1", "#43A047", "#D81B60",
]


@dataclass
class Category:
    """A named, colored set of reactions."""

    name: str
    color: str = "#1A73E8"
    reaction_ids: Set[str] = field(default_factory=set)

    def add(self, reaction_ids) -> None:
        self.reaction_ids.update(reaction_ids)

    def remove(self, reaction_ids) -> None:
        self.reaction_ids.difference_update(reaction_ids)

    def metabolite_ids(self, model: cobra.Model) -> Set[str]:
        """All metabolites touched by this category's reactions."""
        mets: Set[str] = set()
        for rid in self.reaction_ids:
            if model.reactions.has_id(rid):
                mets.update(m.id for m in model.reactions.get_by_id(rid).metabolites)
        return mets

    def to_dict(self) -> dict:
        return {"name": self.name, "color": self.color, "reaction_ids": sorted(self.reaction_ids)}

    @classmethod
    def from_dict(cls, data: dict) -> "Category":
        return cls(
            name=data["name"],
            color=data.get("color", "#1A73E8"),
            reaction_ids=set(data.get("reaction_ids", [])),
        )


class CategoryManager:
    """Holds the categories for a project and assigns colors."""

    def __init__(self) -> None:
        self._categories: Dict[str, Category] = {}
        self._color_index = 0

    def names(self) -> List[str]:
        return list(self._categories.keys())

    def all(self) -> List[Category]:
        return list(self._categories.values())

    def get(self, name: str) -> Category:
        return self._categories[name]

    def has(self, name: str) -> bool:
        return name in self._categories

    def create(self, name: str, color: str | None = None) -> Category:
        if name in self._categories:
            raise ValueError(f"A category named '{name}' already exists.")
        if color is None:
            color = DEFAULT_PALETTE[self._color_index % len(DEFAULT_PALETTE)]
            self._color_index += 1
        cat = Category(name=name, color=color)
        self._categories[name] = cat
        return cat

    def delete(self, name: str) -> None:
        self._categories.pop(name, None)

    def rename(self, old: str, new: str) -> None:
        if new in self._categories:
            raise ValueError(f"A category named '{new}' already exists.")
        cat = self._categories.pop(old)
        cat.name = new
        self._categories[new] = cat

    def category_of_reaction(self, reaction_id: str) -> Category | None:
        for cat in self._categories.values():
            if reaction_id in cat.reaction_ids:
                return cat
        return None

    # persistence -------------------------------------------------------
    def to_list(self) -> list:
        return [c.to_dict() for c in self._categories.values()]

    def load_list(self, data: list) -> None:
        self._categories.clear()
        for entry in data or []:
            cat = Category.from_dict(entry)
            self._categories[cat.name] = cat
        self._color_index = len(self._categories)


def build_subset_model(
    model: cobra.Model,
    reaction_ids,
    *,
    open_boundaries: bool = True,
    boundary_flux: float = 1000.0,
) -> cobra.Model:
    """Extract a standalone sub-model containing only ``reaction_ids``.

    Metabolites shared with reactions outside the subset ("cut" metabolites) get a
    free exchange reaction so the subset can be simulated in isolation, with those
    connections acting as default input/output fluxes.
    """
    keep = {r for r in reaction_ids if model.reactions.has_id(r)}
    if not keep:
        raise ValueError("The category contains no reactions present in the model.")

    sub = model.copy()

    # Identify cut metabolites in the ORIGINAL model: those used by a kept reaction
    # and also by at least one reaction that is NOT kept.
    cut_mets: Set[str] = set()
    for rid in keep:
        for met in model.reactions.get_by_id(rid).metabolites:
            if any(rx.id not in keep for rx in met.reactions):
                cut_mets.add(met.id)

    to_remove = [r for r in sub.reactions if r.id not in keep]
    sub.remove_reactions(to_remove, remove_orphans=True)

    if open_boundaries:
        for mid in cut_mets:
            if not sub.metabolites.has_id(mid):
                continue
            met = sub.metabolites.get_by_id(mid)
            if any(rx.boundary for rx in met.reactions):
                continue
            # Build the exchange reaction manually rather than via add_boundary,
            # whose external-compartment heuristic fails once all the model's
            # original boundary reactions have been removed.
            ex_id = f"EX_subset_{met.id}"
            if sub.reactions.has_id(ex_id):
                continue
            ex = cobra.Reaction(ex_id)
            ex.name = f"Subset exchange for {met.id}"
            ex.lower_bound = -abs(boundary_flux)
            ex.upper_bound = abs(boundary_flux)
            sub.add_reactions([ex])
            ex.add_metabolites({met: -1.0})

    sub.id = f"{model.id or 'model'}_subset"
    return sub
