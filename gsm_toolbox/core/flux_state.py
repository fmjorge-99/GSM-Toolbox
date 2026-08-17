"""Named flux states ("strategies") for the Strategy Explorer.

A *strategy* is a named, solved flux state of the model after a defined set of
edits (knockouts, bound changes, added reactions, objective/medium). Capturing
these lets the graphical engine draw and compare successive rounds of
engineering — colour-coded flux maps, difference maps (Δflux between two
strategies), multi-strategy heatmaps and titre waterfalls.

Pure-Python, no Qt — serialisable so the list of strategies is saved into the
``.gsmtbx`` project (the reproducibility requirement in the proposal).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class FluxState:
    """A named solved flux state (one round of engineering)."""

    name: str
    fluxes: Dict[str, float] = field(default_factory=dict)
    objective_value: float = float("nan")
    method: str = "FBA"
    target: str = ""            # product reaction of interest (for the waterfall)
    notes: str = ""

    def flux(self, reaction_id: str) -> float:
        return float(self.fluxes.get(reaction_id, 0.0))

    def target_flux(self) -> float:
        return abs(self.flux(self.target)) if self.target else float("nan")

    # -- serialisation --------------------------------------------------
    def to_dict(self) -> dict:
        obj = None if self.objective_value != self.objective_value else float(self.objective_value)
        return {
            "name": self.name,
            "fluxes": {k: float(v) for k, v in self.fluxes.items()},
            "objective_value": obj,
            "method": self.method,
            "target": self.target,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "FluxState":
        obj = d.get("objective_value")
        return cls(
            name=d.get("name", "strategy"),
            fluxes={k: float(v) for k, v in (d.get("fluxes") or {}).items()},
            objective_value=float("nan") if obj is None else float(obj),
            method=d.get("method", "FBA"),
            target=d.get("target", ""),
            notes=d.get("notes", ""),
        )


class StrategyManager:
    """An ordered, name-unique collection of :class:`FluxState` strategies."""

    def __init__(self):
        self._states: List[FluxState] = []

    def __len__(self) -> int:
        return len(self._states)

    def __iter__(self):
        return iter(self._states)

    def names(self) -> List[str]:
        return [s.name for s in self._states]

    def get(self, name: str) -> Optional[FluxState]:
        return next((s for s in self._states if s.name == name), None)

    def add(self, state: FluxState) -> FluxState:
        """Add a strategy, making its name unique (append #2, #3… on collision)."""
        base = state.name.strip() or "strategy"
        name, i = base, 2
        existing = set(self.names())
        while name in existing:
            name = f"{base} #{i}"
            i += 1
        state.name = name
        self._states.append(state)
        return state

    def remove(self, name: str) -> None:
        self._states = [s for s in self._states if s.name != name]

    def clear(self) -> None:
        self._states = []

    def difference(self, name_a: str, name_b: str) -> Dict[str, float]:
        """Δflux = state_b − state_a over the union of their reactions.

        Positive values mean the reaction carries *more* flux in ``b`` than in
        ``a`` (drawn red on the diverging map); negative means less (blue)."""
        a, b = self.get(name_a), self.get(name_b)
        if a is None or b is None:
            return {}
        keys = set(a.fluxes) | set(b.fluxes)
        return {k: b.flux(k) - a.flux(k) for k in keys}

    def flux_matrix(self, reaction_ids: List[str]) -> Dict[str, List[float]]:
        """{reaction_id: [flux in each strategy]} for the chosen reactions —
        the data behind the reactions × strategies heatmap / parallel coordinates."""
        return {rid: [s.flux(rid) for s in self._states] for rid in reaction_ids}

    # -- serialisation --------------------------------------------------
    def to_list(self) -> list:
        return [s.to_dict() for s in self._states]

    def load_list(self, data: list) -> None:
        self._states = [FluxState.from_dict(d) for d in (data or [])]
