"""Popup showing a heterologous pathway as a network graph.

The new (heterologous) reactions are highlighted as a category on top of the
host network, together with the native metabolites they connect to. The
embedded network view provides the "steps" expansion: raise the *Steps* control
(or pick a metabolite in *Focus*) to grow the neighbourhood around a chosen node.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import cobra
from PySide6.QtWidgets import QDialog, QLabel, QVBoxLayout

from ..views.network_view import NetworkView


@dataclass
class _PathwayCategory:
    """A minimal stand-in for a Categories entry, understood by NetworkView."""
    name: str
    reaction_ids: List[str] = field(default_factory=list)
    color: str = "#E86A5C"


class PathwayGraphDialog(QDialog):
    """Show a pathway (host + new reactions) in an interactive network view."""

    def __init__(self, parent, model: cobra.Model, reaction_ids: List[str],
                 target: str = "", category_name: str = "New pathway",
                 added_metabolites: List[str] | None = None):
        super().__init__(parent)
        self.setWindowTitle(f"Pathway map — {target}" if target else "Pathway map")
        self.resize(1000, 700)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

        layout = QVBoxLayout(self)
        caption = QLabel(
            "Newly designed reactions are highlighted. "
            "<span style='color:#34A853'><b>Green</b></span> metabolites are not native "
            "to the model; the legend below maps the remaining colours to "
            "compartments. Raise <b>Steps</b> to expand the neighbourhood.")
        caption.setWordWrap(True)
        layout.addWidget(caption)

        self.view = NetworkView()
        self.view.set_render_busy(False)   # embedded view: no nested busy popup (#B8)
        layout.addWidget(self.view, 1)

        # Only reactions actually present in the model can be highlighted.
        present = [rid for rid in reaction_ids if model.reactions.has_id(rid)]
        cat = _PathwayCategory(name=category_name, reaction_ids=present)
        self.view.set_model(model)
        # Green = came from the database (created de novo), blue = native precursor.
        # Without this every metabolite rendered the same blue, so there was no way to
        # tell what the host already makes from what the pathway introduces.
        self.view.set_added_metabolites(added_metabolites or [])
        self.view.set_categories([cat])
        # Show just the pathway + directly-connected metabolites first (Steps 0).
        self.view.radius_spin.setValue(0)
        self.view.focus_category(category_name)
