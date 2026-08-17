"""Plot Gallery — pick a publication-grade figure and export it.

A light host around the shared PlotCanvas: it takes a set of named render
callbacks (assembled by the main window from the current analysis result, the
last FBA fluxes and the saved strategies), shows them in a dropdown, and renders
the chosen one with figure (SVG/PDF/PNG) + data (CSV) export.
"""

from __future__ import annotations

from typing import Callable, List, Tuple

from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, QDialog

from ..viz.plot_canvas import PlotCanvas
from ..widgets.dialog_util import clamp_to_screen


class PlotGalleryDialog(QDialog):
    def __init__(self, parent, entries: List[Tuple[str, Callable[[PlotCanvas], None]]]):
        super().__init__(parent)
        self.setWindowTitle("Plot Gallery")
        self.resize(880, 620)
        self._entries = entries

        self.combo = QComboBox()
        self.combo.addItems([label for label, _ in entries])
        self.combo.currentIndexChanged.connect(self._render)
        top = QHBoxLayout()
        top.addWidget(QLabel("Figure:"))
        top.addWidget(self.combo, 1)

        self.canvas = PlotCanvas()
        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.canvas, 1)

        clamp_to_screen(self)
        if entries:
            self._render()

    def _render(self) -> None:
        i = self.combo.currentIndex()
        if 0 <= i < len(self._entries):
            self._entries[i][1](self.canvas)
