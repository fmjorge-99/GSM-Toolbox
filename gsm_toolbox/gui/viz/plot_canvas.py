"""A reusable matplotlib canvas with publication export.

Wraps a Figure + navigation toolbar and adds "Export figure…" (SVG/PDF/PNG) and
"Export data (CSV)" — so every graphical-engine plot exports both a vector figure
and the values behind it, per the proposal's reproducibility requirement.
"""

from __future__ import annotations

from typing import Optional

import matplotlib

matplotlib.use("QtAgg")

import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PySide6.QtWidgets import QHBoxLayout, QMessageBox, QPushButton, QVBoxLayout, QWidget

from . import theme
from ..widgets.dialog_util import choose_save_path


class PlotCanvas(QWidget):
    def __init__(self, parent=None, *, figsize=(6.4, 4.4)):
        super().__init__(parent)
        theme.apply_style()
        self.figure = Figure(figsize=figsize, tight_layout=True)
        self.figure.set_facecolor("white")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self._data: Optional[pd.DataFrame] = None
        self._title = "figure"

        self.export_fig_btn = QPushButton("Export figure…")
        self.export_fig_btn.setToolTip("Save the figure as vector SVG/PDF (or PNG).")
        self.export_fig_btn.clicked.connect(self._export_figure)
        self.export_csv_btn = QPushButton("Export data (CSV)")
        self.export_csv_btn.setToolTip("Save the values behind this figure as CSV.")
        self.export_csv_btn.clicked.connect(self._export_csv)
        btn_row = QHBoxLayout()
        btn_row.addWidget(self.toolbar, 1)
        btn_row.addWidget(self.export_fig_btn)
        btn_row.addWidget(self.export_csv_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(btn_row)
        layout.addWidget(self.canvas, 1)
        self._update_enabled()

    def render(self, builder, *args, title: str = "figure", **kwargs) -> None:
        """Call a plots.* builder(fig, *args) that draws into our Figure and
        returns the underlying DataFrame; then refresh and remember it for export."""
        self._title = title
        try:
            self._data = builder(self.figure, *args, **kwargs)
        except Exception as exc:  # noqa: BLE001 - never let a plot crash the app
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, f"Could not draw this plot:\n{exc}", ha="center",
                    va="center", color="#D93025", fontsize=9, wrap=True)
            self._data = None
        self.canvas.draw()
        self._update_enabled()

    def _update_enabled(self) -> None:
        self.export_csv_btn.setEnabled(self._data is not None and not self._data.empty)

    def _export_figure(self) -> None:
        path = choose_save_path(self, "Export figure", f"{self._title}.svg",
                                "Vector SVG (*.svg);;PDF (*.pdf);;PNG image (*.png)")
        if not path:
            return
        if not path.lower().endswith((".svg", ".pdf", ".png")):
            path += ".svg"
        try:
            self.figure.savefig(path, dpi=200, bbox_inches="tight")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Could not export figure", str(exc))

    def _export_csv(self) -> None:
        if self._data is None or self._data.empty:
            return
        path = choose_save_path(self, "Export data", f"{self._title}.csv", "CSV (*.csv)")
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        try:
            self._data.to_csv(path, index=False)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Could not export data", str(exc))
