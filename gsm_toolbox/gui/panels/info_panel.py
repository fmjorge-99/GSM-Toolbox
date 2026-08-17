"""Right-dock information panel: a 'Model' tab (general info) and a 'Selection' tab.

Combines the model summary dashboard and the per-item inspector into one tabbed
widget, and auto-switches to the Selection tab when the user picks something.
"""

from __future__ import annotations

from typing import Optional

import cobra
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QTabWidget

from ..views.summary_view import SummaryView
from .inspector import Inspector


class InfoPanel(QTabWidget):
    enzymes_requested = Signal(list, str)  # [ec, …], reaction id

    def __init__(self):
        super().__init__()
        self.summary = SummaryView()
        self.inspector = Inspector()
        self.inspector.enzymes_requested.connect(self.enzymes_requested)
        self.addTab(self.summary, "Model")
        self.addTab(self.inspector, "Selection")
        # Allow the panel to be made narrow / short so it never crowds the center.
        self.setMinimumWidth(150)
        self.setMinimumHeight(90)

    # model-level info --------------------------------------------------
    def update_model(self, model: Optional[cobra.Model], diff: dict | None = None,
                     last_growth: Optional[float] = None) -> None:
        self.summary.update_summary(model, diff, last_growth)

    # selection info (auto-focus the Selection tab) ---------------------
    def show_reaction(self, rxn: cobra.Reaction, flux=None) -> None:
        self.inspector.show_reaction(rxn, flux)
        self.setCurrentWidget(self.inspector)

    def show_metabolite(self, met: cobra.Metabolite) -> None:
        self.inspector.show_metabolite(met)
        self.setCurrentWidget(self.inspector)

    def show_gene(self, gene: cobra.Gene) -> None:
        self.inspector.show_gene(gene)
        self.setCurrentWidget(self.inspector)

    def clear_selection(self) -> None:
        self.inspector.show_placeholder()
