"""Embedded matplotlib canvas for analysis plots (envelopes, robustness, PhPP)."""

from __future__ import annotations

import matplotlib

matplotlib.use("QtAgg")  # ensure the Qt backend before importing pyplot machinery

import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PySide6.QtWidgets import QVBoxLayout, QWidget

from .. import style


class PlotView(QWidget):
    def __init__(self):
        super().__init__()
        self.figure = Figure(figsize=(5, 4), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.title = "plot"
        layout = QVBoxLayout(self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, 1)

    def _fresh_axes(self, title: str = "plot"):
        self.title = title
        self.figure.clear()
        return self.figure.add_subplot(111)

    def plot_production_envelope(self, df: pd.DataFrame, target: str) -> None:
        ax = self._fresh_axes(f"production_envelope_{target}")
        # cobra's production_envelope has flux_minimum/flux_maximum vs the target flux.
        x = df[target] if target in df.columns else df.iloc[:, -1]
        ax.fill_between(x, df["flux_minimum"], df["flux_maximum"],
                        color=style.ACCENT, alpha=0.25, label="Feasible range")
        ax.plot(x, df["flux_maximum"], color=style.ACCENT, lw=2, label="Max growth")
        ax.plot(x, df["flux_minimum"], color=style.FLUX_REVERSE, lw=1.5, label="Min growth")
        ax.set_xlabel(f"{target} flux")
        ax.set_ylabel("Growth rate")
        ax.set_title("Production envelope")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        self.canvas.draw()

    def plot_robustness(self, df: pd.DataFrame, control: str) -> None:
        ax = self._fresh_axes(f"robustness_{control}")
        clean = df.dropna(subset=["objective"])  # drop infeasible fixed-flux points
        ax.plot(clean["control_flux"], clean["objective"], color=style.ACCENT,
                lw=2, marker="o", ms=3)
        ax.set_xlabel(f"{control} flux (fixed)")
        ax.set_ylabel("Objective (growth)")
        ax.set_title("Robustness analysis")
        ax.grid(alpha=0.3)
        if len(clean) <= 1:
            ax.text(0.5, 0.5, "Only one feasible point in this range —\n"
                    "widen the scan range or disable auto-range.",
                    ha="center", va="center", transform=ax.transAxes, color=style.TEXT_MUTED)
        self.canvas.draw()

    def plot_phase_plane(self, df: pd.DataFrame, x_name: str, y_name: str) -> None:
        ax = self._fresh_axes(f"phase_plane_{x_name}_{y_name}")
        pivot = df.pivot(index="y", columns="x", values="objective")
        im = ax.imshow(
            pivot.values, origin="lower", aspect="auto", cmap="viridis",
            extent=[pivot.columns.min(), pivot.columns.max(),
                    pivot.index.min(), pivot.index.max()],
        )
        ax.set_xlabel(f"{x_name} flux")
        ax.set_ylabel(f"{y_name} flux")
        ax.set_title("Phenotypic phase plane (objective)")
        self.figure.colorbar(im, ax=ax, label="Objective")
        self.canvas.draw()

    def plot_line(self, df: pd.DataFrame, x: str, y: str, *, title: str,
                  xlabel: str, ylabel: str) -> None:
        ax = self._fresh_axes(title)
        ax.plot(df[x], df[y], color=style.ACCENT, lw=2, marker="o", ms=3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.3)
        self.canvas.draw()

    def has_content(self) -> bool:
        return len(self.figure.axes) > 0

    def save(self, path: str) -> None:
        """Save the current figure to ``path`` (format inferred from extension)."""
        self.figure.savefig(path, dpi=200, bbox_inches="tight")

    def clear(self) -> None:
        self.figure.clear()
        self.canvas.draw()
