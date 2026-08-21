"""Show flux-analysis results as a network map.

Reactions carrying flux are drawn as a focused sub-network; each connector is
coloured and thickened by the flux magnitude/direction (as in published FBA
figures) and each reaction node is labelled with its flux value — and, for FVA,
its [lower, upper] bounds. Built on the standard NetworkView renderer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QLabel, QVBoxLayout

from ..views.network_view import NetworkView


@dataclass
class _FluxCategory:
    name: str
    reaction_ids: List[str] = field(default_factory=list)
    color: str = "#1a73e8"


class FluxMapDialog(QDialog):
    def __init__(self, parent, model: cobra.Model, reaction_ids: List[str],
                 fluxes: Dict[str, float], labels: Dict[str, str]):
        super().__init__(parent)
        self.setWindowTitle("Flux network")
        self.resize(1000, 720)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

        layout = QVBoxLayout(self)
        caption = QLabel(
            "Connector width and colour scale with flux magnitude and direction; each "
            "reaction is labelled with its value. Use <b>Steps</b> to expand.")
        caption.setWordWrap(True)
        layout.addWidget(caption)

        self.view = NetworkView()
        self.view.set_render_busy(False)   # embedded view: no nested busy popup (#B8)
        layout.addWidget(self.view, 1)

        cat = _FluxCategory(name="Flux-carrying reactions", reaction_ids=reaction_ids)
        self.view.set_model(model)
        self.view.set_categories([cat])
        self.view.set_fluxes(fluxes)
        self.view.set_flux_values(labels)
        self.view.currency_check.setChecked(True)
        self.view.radius_spin.setValue(0)
        self.view.focus_category("Flux-carrying reactions")
