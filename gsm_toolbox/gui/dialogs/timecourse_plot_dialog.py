"""Plot a time course: pick the axes, look at it, export it.

A dynamic run produces a wide table — time, biomass, growth rate, one column per medium
component and one per reported flux. The interesting relationship is rarely the whole
table: it is biomass against time, or a flux against the nutrient that is running out.
Reading that off a grid of numbers is exactly the step where a phase change gets missed.

The X axis is free rather than fixed to time on purpose. Plotting a flux against the
remaining nitrate answers "what does the cell do as nitrogen runs out" directly, without
the reader having to hold two columns in their head at once.

Series on different scales (biomass ~0.1 gDW L⁻¹, a flux ~20 mmol gDW⁻¹ h⁻¹) are put on a
second axis rather than being silently flattened into one — a shared axis would render the
smaller series as a line along zero and invite the conclusion that nothing happened.
"""
from __future__ import annotations

from typing import List, Optional

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel,
    QListWidget, QListWidgetItem, QMessageBox, QPushButton, QVBoxLayout, QWidget)

from .. import style
from ..views.plot_view import PlotView

#: Columns that describe the run rather than measure it — never offered as a series.
_NON_NUMERIC = {"rules_fired"}

#: Beyond this many series, name them in the legend rather than on the axis.
_MAX_NAMED_ON_AXIS = 3


def _axis_label(columns) -> str:
    if not columns:
        return ""
    if len(columns) <= _MAX_NAMED_ON_AXIS:
        return ", ".join(columns)
    return f"{len(columns)} series (see legend)"


class TimecoursePlotDialog(QDialog):
    """Choose X and Y series from a time-course table and render them."""

    def __init__(self, parent, frame: pd.DataFrame, title: str = "Time course",
                 runs: Optional[List[tuple]] = None):
        """``runs`` is [(name, frame), …]; several may be overlaid on one plot.

        Overlaying is the whole reason runs are stored separately. A regulated run and an
        unregulated one differ by a few percent in biomass and completely in when a
        transition happens — a difference that is obvious on one pair of axes and nearly
        invisible across two windows.
        """
        super().__init__(parent)
        self.setWindowTitle(f"Plot — {title}")
        self.resize(1080, 660)
        from ..widgets.dialog_util import clamp_to_screen

        self._runs: List[tuple] = list(runs) if runs else [(title, frame)]
        self._frame = frame if frame is not None else self._runs[0][1]
        self._title = title

        outer = QVBoxLayout(self)
        hint = QLabel(
            "Choose what goes on each axis. Plotting a flux against a nutrient "
            "concentration, rather than against time, shows the response to the nutrient "
            "directly. Select more than one run to overlay them.")
        hint.setWordWrap(True)
        outer.addWidget(hint)

        body = QHBoxLayout()
        body.addWidget(self._controls(), 0)
        self.plot = PlotView()
        body.addWidget(self.plot, 1)
        outer.addLayout(body, 1)

        buttons = QHBoxLayout()
        self.draw_btn = QPushButton("Draw")
        self.draw_btn.setObjectName("primary")
        self.draw_btn.clicked.connect(self.draw)
        buttons.addWidget(self.draw_btn)
        export = QPushButton("Export image…")
        export.clicked.connect(self._export_image)
        buttons.addWidget(export)
        export_data = QPushButton("Export data…")
        export_data.setToolTip("Save the plotted columns as CSV.")
        export_data.clicked.connect(self._export_data)
        buttons.addWidget(export_data)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        outer.addLayout(buttons)

        self.draw()
        clamp_to_screen(self)

    # -- controls -------------------------------------------------------------------
    def _selected_runs(self) -> List[tuple]:
        chosen = {i.text() for i in self.run_list.selectedItems()} \
            if getattr(self, "run_list", None) else set()
        picked = [(n, f) for n, f in self._runs if n in chosen]
        return picked or self._runs[:1]

    def _numeric_columns(self) -> List[str]:
        """Columns offered as axes: the union across the selected runs.

        A union rather than an intersection, because two runs of the same kind can
        legitimately report different fluxes and hiding a column that exists in one of
        them would silently narrow what can be asked. A run missing a column is simply
        not drawn for it, and the note says so.
        """
        out: List[str] = []
        for _name, frame in self._selected_runs():
            if frame is None:
                continue
            for column in frame.columns:
                if column in _NON_NUMERIC or str(column) in out:
                    continue
                if pd.api.types.is_numeric_dtype(frame[column]):
                    out.append(str(column))
        return out

    def _controls(self) -> QWidget:
        panel = QWidget()
        panel.setMaximumWidth(340)
        v = QVBoxLayout(panel)
        v.setContentsMargins(0, 0, 0, 0)

        # Runs first: what is selected here decides which columns the axes can offer.
        self.run_list = QListWidget()
        self.run_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        for name, _frame in self._runs:
            self.run_list.addItem(QListWidgetItem(name))
        if self._runs:
            self.run_list.item(0).setSelected(True)
        if len(self._runs) > 1:
            v.addWidget(QLabel("<b>Runs</b> (select several to overlay)"))
            self.run_list.setMaximumHeight(110)
            v.addWidget(self.run_list)
            self.run_list.itemSelectionChanged.connect(self._sync_run_selection)
        else:
            self.run_list.setVisible(False)

        columns = self._numeric_columns()

        v.addWidget(QLabel("<b>X axis</b>"))
        self.x_axis = QComboBox()
        for column in columns:
            self.x_axis.addItem(column)
        if "time_h" in columns:
            self.x_axis.setCurrentText("time_h")
        v.addWidget(self.x_axis)

        v.addWidget(QLabel("<b>Y axis</b> (one or more)"))
        self.y_axis = QListWidget()
        self.y_axis.setSelectionMode(QAbstractItemView.ExtendedSelection)
        for column in columns:
            if column == self.x_axis.currentText():
                continue
            self.y_axis.addItem(QListWidgetItem(column))
        v.addWidget(self.y_axis, 1)
        # Biomass is what a batch run is usually about; preselecting it means the dialog
        # opens on a useful plot rather than an empty frame.
        for preferred in ("biomass_gDW_L", "growth_rate_per_h"):
            items = self.y_axis.findItems(preferred, Qt.MatchExactly)
            if items:
                items[0].setSelected(True)
                break
        self.x_axis.currentTextChanged.connect(self._sync_y_options)

        self.second_axis = QCheckBox("Separate axis for differing scales")
        self.second_axis.setChecked(True)
        self.second_axis.setToolTip(
            "Series whose magnitudes differ by more than 20× get their own axis, so a "
            "small series is not drawn as a flat line at zero.")
        v.addWidget(self.second_axis)

        self.markers = QCheckBox("Show points")
        self.markers.setChecked(True)
        self.markers.setToolTip(
            "Points show where the solver actually evaluated — useful for judging "
            "whether the integration step is fine enough.")
        v.addWidget(self.markers)

        self.note = QLabel("")
        self.note.setWordWrap(True)
        self.note.setStyleSheet(f"color:{style.TEXT_MUTED};")
        v.addWidget(self.note)
        return panel

    def _sync_y_options(self) -> None:
        chosen = {i.text() for i in self.y_axis.selectedItems()}
        self.y_axis.clear()
        for column in self._numeric_columns():
            if column == self.x_axis.currentText():
                continue
            item = QListWidgetItem(column)
            self.y_axis.addItem(item)
            if column in chosen:
                item.setSelected(True)

    def _sync_run_selection(self) -> None:
        """Changing which runs are shown can change which columns exist."""
        columns = self._numeric_columns()
        current_x = self.x_axis.currentText()
        self.x_axis.blockSignals(True)
        self.x_axis.clear()
        for column in columns:
            self.x_axis.addItem(column)
        if current_x in columns:
            self.x_axis.setCurrentText(current_x)
        elif "time_h" in columns:
            self.x_axis.setCurrentText("time_h")
        self.x_axis.blockSignals(False)
        self._sync_y_options()
        self.draw()

    # -- drawing --------------------------------------------------------------------
    def _selected_y(self) -> List[str]:
        return [i.text() for i in self.y_axis.selectedItems()]

    def draw(self) -> None:
        x = self.x_axis.currentText()
        ys = self._selected_y()
        axes = self.plot._fresh_axes(f"timecourse_{self._title}")
        if not x or not ys:
            axes.text(0.5, 0.5, "Choose an X axis and at least one Y series.",
                      ha="center", va="center", transform=axes.transAxes, fontsize=10)
            axes.set_axis_off()
            self.plot.canvas.draw()
            self.note.setText("")
            return

        runs = self._selected_runs()
        marker = "o" if self.markers.isChecked() else None
        multi = len(runs) > 1

        # Which run/column pairs actually exist. A column missing from one run is
        # skipped for that run and reported, never drawn as a line at zero.
        pairs, missing = [], []
        for name, frame in runs:
            if frame is None or x not in frame.columns:
                missing.append(f"{name} (no “{x}”)")
                continue
            ordered = frame.sort_values(by=x)
            for column in ys:
                if column in ordered.columns:
                    pairs.append((name, column, ordered))
                else:
                    missing.append(f"{name}·{column}")
        if not pairs:
            axes.text(0.5, 0.5, "None of the selected runs has these columns.",
                      ha="center", va="center", transform=axes.transAxes, fontsize=10)
            axes.set_axis_off()
            self.plot.canvas.draw()
            self.note.setText("; ".join(missing[:6]))
            return

        # Split by magnitude so a small series is not flattened against a large one.
        # Scale is judged per column across all runs, so the same quantity keeps the
        # same axis in every run and the overlay stays comparable.
        scales: dict = {}
        for _name, column, frame in pairs:
            peak = float(pd.to_numeric(frame[column], errors="coerce").abs().max() or 0.0)
            scales[column] = max(scales.get(column, 0.0), peak)
        biggest = max(scales.values()) if scales else 0.0
        secondary = set()
        if self.second_axis.isChecked() and biggest > 0:
            secondary = {c for c, v in scales.items() if v > 0 and biggest / v > 20}

        palette = style.SERIES_COLORS
        column_colour = {c: palette[i % len(palette)] for i, c in enumerate(ys)}
        column_marker = {c: style.SERIES_MARKERS[i % len(style.SERIES_MARKERS)]
                         for i, c in enumerate(ys)}
        if multi:
            # Colour is the quantity, dash is the run: the eye reads "same measurement"
            # from the colour and "which run" from the pattern, which is the comparison
            # an overlay exists to make.
            line_style = {(name, column): style.SERIES_DASHES[i % len(style.SERIES_DASHES)]
                          for i, (name, _f) in enumerate(runs) for column in ys}
        else:
            # One run: the dash is free, so spend it on the series. Twelve colours are
            # mutually distinguishable, but a plot may follow more reactions than that,
            # and past twelve a repeat of colour alone is exactly the complaint this
            # palette was chosen to fix.
            line_style = {(runs[0][0], column):
                          style.SERIES_DASHES[(i // len(palette)) % len(style.SERIES_DASHES)]
                          for i, column in enumerate(ys)}

        twin = None
        handles, labels = [], []
        for name, column, frame in pairs:
            target = axes
            if column in secondary:
                twin = twin or axes.twinx()
                target = twin
            label = f"{name} · {column}" if multi else column
            line, = target.plot(
                frame[x], frame[column], lw=2, ms=5,
                marker=column_marker.get(column) if marker else None,
                markevery=max(1, len(frame) // 12),
                color=column_colour.get(column, style.SERIES_COLORS[0]),
                ls=line_style.get((name, column), "-"),
                label=label)
            handles.append(line)
            labels.append(label)

        primary_cols = [c for c in ys if c not in secondary]
        axes.set_xlabel(x)
        # Past a few series the joined names are longer than the axis and spill off the
        # figure; the legend already carries identity, so the axis just says what it is.
        axes.set_ylabel(_axis_label(primary_cols))
        axes.grid(alpha=0.3)
        notes = []
        if twin is not None:
            twin.set_ylabel(_axis_label([c for c in ys if c in secondary]))
            notes.append("Right-hand axis carries the series whose scale differs too "
                         "much to share one.")
        if multi:
            notes.append("Colour is the quantity, line style is the run.")
        if missing:
            notes.append("Not present: " + ", ".join(missing[:6])
                         + (" …" if len(missing) > 6 else "") + ".")
        self.note.setText(" ".join(notes))

        if handles:
            axes.legend(handles, labels, fontsize=8, loc="best")
        axes.set_title(" · ".join(n for n, _f in runs) if multi else self._title)
        self.plot.figure.tight_layout()
        self.plot.canvas.draw()

    # -- export ---------------------------------------------------------------------
    def _export_image(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export plot", f"{self._title}.png",
            "PNG image (*.png);;SVG vector (*.svg);;PDF (*.pdf)")
        if not path:
            return
        try:
            self.plot.figure.savefig(path, dpi=200, bbox_inches="tight")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Exported", f"Plot written to\n{path}")

    def _export_data(self) -> None:
        """Export exactly what is plotted — every selected run, tagged by name."""
        wanted = [self.x_axis.currentText()] + self._selected_y()
        pieces = []
        for name, frame in self._selected_runs():
            if frame is None:
                continue
            columns = [c for c in wanted if c in frame.columns]
            if not columns:
                continue
            piece = frame[columns].copy()
            piece.insert(0, "run", name)
            pieces.append(piece)
        if not pieces:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export plotted data", f"{self._title}.csv", "CSV (*.csv)")
        if not path:
            return
        try:
            # Long form, one row per run per point: it keeps runs with different columns
            # or different lengths in one file without inventing blanks to align them.
            pd.concat(pieces, ignore_index=True).to_csv(path, index=False)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Exported", f"Data written to\n{path}")
