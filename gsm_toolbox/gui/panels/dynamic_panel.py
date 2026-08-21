"""Dynamic Analysis — condition scanning and time-course simulation.

The Analysis tab answers steady-state questions: what can the cell do under one fixed set
of conditions. Most regulation, though, is about what happens as conditions *change* —
nitrogen running out over a batch, inorganic carbon drawn down by a dense culture, light
switching between day and night. Those need a different shape of experiment, and mixing
them into the steady-state tab would blur a distinction that matters for how the results
should be read.

Two modes:

* **Condition scan** — sweep one variable across a range, solve at each point, and show
  where the network rewires. The interesting part of a scan is its breakpoints, so those
  are called out rather than left for the reader to spot in a table.
* **Time course** — dynamic FBA over a batch: the chosen medium components deplete, the
  regulatory state is re-evaluated from the *current* concentrations, and the culture
  switches phase when it crosses a threshold.

Both work from **the model's own medium** rather than a fixed list of nutrients. A
hard-coded "nitrogen / carbon / iron" menu cannot express the experiment most users
actually want to run, and it silently misrepresents any model whose medium differs. What
keeps that generality biologically meaningful is the sensor indirection: each component
still maps onto a named sensor (nitrogen, inorganic carbon…), so a rule written for
nitrogen limitation keeps working when the culture is switched from nitrate to ammonium.
"""
from __future__ import annotations

import re
from collections import OrderedDict
from typing import Dict, List, Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QGroupBox, QFormLayout,
    QHBoxLayout, QHeaderView, QInputDialog, QLabel, QMessageBox, QPushButton,
    QScrollArea, QSizePolicy, QSpinBox, QTabBar, QTabWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

from ...core import regulation as reg
from .. import style

#: Sensors offered as scan variables even when no exchange supplies them directly.
#: These are the organism-level quantities a rule set is written against.
BASE_VARIABLES = [
    ("light_uE", "Light (µmol photons m⁻² s⁻¹)", 0.0, 3000.0),
    ("inorganic_c_mM", "Inorganic carbon (mM)", 0.0, 100.0),
    ("nitrogen_mM", "Nitrogen, any source (mM)", 0.0, 100.0),
    ("iron_uM", "Iron (µM)", 0.0, 100.0),
]

#: Starting concentrations offered for well-known components, in mM. Anything not listed
#: starts at 10 mM, which is a placeholder the user is expected to set — the table shows
#: it plainly rather than hiding a guess inside the solver.
DEFAULT_CONCENTRATION = {
    "no3": 17.6,    # BG-11
    "nh4": 5.0,
    "hco3": 50.0,
    "co2": 0.013,
    "pi": 0.18,
    "so4": 0.3,
    "glc__D": 5.0,
    "urea": 5.0,
}

_MODE_DEPLETE = "Depletes"
_MODE_BUFFERED = "Held constant"


def _elide(combo: QComboBox, characters: int = 18) -> None:
    """Stop a combo sizing itself to its longest entry.

    Left alone, a box holding "Light (µmol photons m⁻² s⁻¹)" claims ~370 px of minimum
    width. A few of those is all it takes to push the main window's minimum past a
    1366-px laptop screen, which Qt resolves by dropping the window out of Maximized —
    the recurring layout bug this project has a guard test for.
    """
    combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
    combo.setMinimumContentsLength(characters)


class DynamicAnalysisPanel(QWidget):
    """Condition scanning and time-course simulation, side by side with their results."""

    scan_requested = Signal(dict)
    timecourse_requested = Signal(dict)
    cancel_requested = Signal()
    #: Open the rule-set manager/editor without leaving this tab.
    regulation_requested = Signal()
    #: Plot the runs the user has stored, by name.
    plot_requested = Signal(list)
    #: Show one stored run's table.
    table_requested = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._model = None
        self._scan_targets: List[str] = []
        self._tc_targets: List[str] = []
        self._last_kind = ""
        #: name -> {frame, title, commentary, warning, kind}. Ordered as the tabs are.
        self._runs: "OrderedDict[str, dict]" = OrderedDict()
        self._counters = {"scan": 0, "timecourse": 0}

        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(6)

        # The settings take the whole panel. Results no longer share the space with
        # them: a run is stored under a tab and opened in its own window, which is what
        # frees the parameter rows to be legible instead of compressed into a strip.
        #
        # Note there is still no splitter and no manual sizing. An earlier version
        # recomputed a height inside ``resizeEvent`` and applied ``setFixedHeight``;
        # setting the height triggered a relayout, the relayout delivered another resize
        # event, and on a maximized window the two fed each other until the application
        # stopped responding. Qt resolves this layout in one pass with no code of ours
        # running during it.
        self.controls = QTabWidget()
        self.controls.addTab(self._build_scan_tab(), "Condition scan")
        self.controls.addTab(self._build_timecourse_tab(), "Time course")
        self.controls.currentChanged.connect(self._sync_run_button)
        outer.addWidget(self.controls, 1)

        outer.addLayout(self._action_row(), 0)
        outer.addWidget(self._run_tabs(), 0)
        self._sync_run_button()

    @staticmethod
    def _scrollable(content: QWidget) -> QScrollArea:
        area = QScrollArea()
        area.setWidget(content)
        area.setWidgetResizable(True)
        area.setFrameShape(QScrollArea.NoFrame)
        # Horizontal scrolling would mean the layout is too wide, which is the window-width
        # bug in another guise; the controls are sized to fit, so only vertical is offered.
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        return area

    def _page(self, form: QWidget, hint: str) -> QWidget:
        """A tab page: scrollable settings with a one-line note beneath them."""
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 0, 4)
        layout.setSpacing(4)
        layout.addWidget(self._scrollable(form), 1)
        note = self._muted(hint)
        note.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        layout.addWidget(note, 0)
        return page

    # -- actions -------------------------------------------------------------------
    def _action_row(self) -> QHBoxLayout:
        """Run, and everything you do with what a run produced."""
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)

        self.run_btn = QPushButton("Run")
        self.run_btn.setObjectName("primary")
        self.run_btn.clicked.connect(self._emit_run)
        row.addWidget(self.run_btn, 0)

        self.plot_btn = QPushButton("Plot Results")
        self.plot_btn.setToolTip("Choose the axes and draw. Several stored runs can be "
                                 "overlaid on one plot.")
        self.plot_btn.clicked.connect(self._emit_plot)
        row.addWidget(self.plot_btn, 0)

        self.table_btn = QPushButton("Display Table")
        self.table_btn.setToolTip("Open the selected run's results table.")
        self.table_btn.clicked.connect(self._emit_table)
        row.addWidget(self.table_btn, 0)

        row.addSpacing(12)
        self.use_regulation = QCheckBox("Regulation")
        self.use_regulation.setToolTip(
            "Apply the active rule set to the run. Off runs the plain model.")
        row.addWidget(self.use_regulation, 0)

        self.rules_btn = QPushButton("Regulation settings…")
        self.rules_btn.clicked.connect(self.regulation_requested.emit)
        row.addWidget(self.rules_btn, 0)

        self.rules_status = QLabel("")
        self.rules_status.setStyleSheet(f"color:{style.TEXT_MUTED};")
        self.rules_status.setTextFormat(Qt.RichText)
        # Shrinkable, with the full text in the tooltip: a status line that refuses to
        # give way would push the window's minimum width past the screen.
        self.rules_status.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        row.addWidget(self.rules_status, 1)

        self._sync_result_buttons()
        return row

    # -- stored runs ---------------------------------------------------------------
    def _run_tabs(self) -> QWidget:
        """One tab per completed run: closable, renamable, and the plot's selection.

        Runs accumulate rather than replacing one another because the comparison between
        two runs — regulated against unregulated, one medium against another — is the
        thing being looked for, and it cannot be seen if the second overwrites the first.
        """
        self.run_tabs = QTabBar()
        self.run_tabs.setTabsClosable(True)
        self.run_tabs.setMovable(True)
        self.run_tabs.setExpanding(False)
        self.run_tabs.setDrawBase(False)
        self.run_tabs.setUsesScrollButtons(True)
        self.run_tabs.tabCloseRequested.connect(self._close_run)
        self.run_tabs.tabBarDoubleClicked.connect(self._rename_run)
        self.run_tabs.currentChanged.connect(lambda _i: self._sync_result_buttons())

        self.run_hint = self._muted(
            "Completed runs appear here. Double-click a tab to rename it.")
        self.run_hint.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)

        holder = QWidget()
        layout = QHBoxLayout(holder)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(self.run_tabs, 0)
        layout.addWidget(self.run_hint, 1)
        return holder

    # -- controls ------------------------------------------------------------------
    def _build_scan_tab(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)

        box = QGroupBox("Sweep one condition")
        form = QFormLayout(box)
        self.scan_variable = QComboBox()
        _elide(self.scan_variable)
        self.scan_variable.setToolTip(
            "A named sensor, or any exchange in the model. Exchanges are converted to "
            "an uptake bound; sensors are read directly by the rule set.")
        self.scan_variable.currentIndexChanged.connect(self._sync_scan_range)
        form.addRow("Variable:", self.scan_variable)

        row = QHBoxLayout()
        self.scan_from = QDoubleSpinBox()
        self.scan_to = QDoubleSpinBox()
        for w in (self.scan_from, self.scan_to):
            w.setDecimals(3)
            w.setRange(0.0, 100000.0)
            w.setMaximumWidth(110)
        row.addWidget(self.scan_from)
        row.addWidget(QLabel("→"))
        row.addWidget(self.scan_to)
        row.addStretch(1)
        form.addRow("Range:", row)

        self.scan_points = QSpinBox()
        self.scan_points.setRange(2, 60)
        self.scan_points.setValue(10)
        self.scan_points.setMaximumWidth(70)
        form.addRow("Steps:", self.scan_points)

        self.scan_log = QComboBox()
        self.scan_log.addItem("Linear spacing", False)
        self.scan_log.addItem("Logarithmic spacing", True)
        _elide(self.scan_log)
        self.scan_log.setToolTip(
            "Logarithmic spacing is usually right for a nutrient concentration, where "
            "the interesting behaviour is near zero.")
        form.addRow("Spacing:", self.scan_log)

        target_row = QHBoxLayout()
        self.scan_targets_btn = QPushButton("Choose reactions…")
        self.scan_targets_btn.setToolTip(
            "Pick which fluxes are recorded at each point, besides growth.")
        self.scan_targets_btn.clicked.connect(lambda: self._pick_targets("scan"))
        target_row.addWidget(self.scan_targets_btn)
        self.scan_targets_label = QLabel("growth only")
        self.scan_targets_label.setStyleSheet(f"color:{style.TEXT_MUTED};")
        target_row.addWidget(self.scan_targets_label, 1)
        form.addRow("Report:", target_row)
        self._fix_height(self.scan_variable, self.scan_from, self.scan_to,
                         self.scan_points, self.scan_log, self.scan_targets_btn)
        v.addWidget(box)

        return self._page(
            page,
            "Breakpoints — where the response changes sharply — are reported with the "
            "results.")

    def _build_timecourse_tab(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)

        box = QGroupBox("Medium components to follow")
        bv = QVBoxLayout(box)
        self.substrates = QTableWidget(0, 5)
        self.substrates.setHorizontalHeaderLabels(
            ["Component", "Exchange", "Initial (mM)", "Max uptake", "Mode"])
        self.substrates.verticalHeader().setVisible(False)
        self.substrates.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.substrates.setMinimumHeight(120)
        header = self.substrates.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        for c in (1, 2, 3, 4):
            header.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        bv.addWidget(self.substrates)

        srow = QHBoxLayout()
        self.add_substrate_btn = QPushButton("Add component…")
        self.add_substrate_btn.setToolTip(
            "Any exchange the model has. Components in the medium are listed first.")
        self.add_substrate_btn.clicked.connect(self._pick_substrates)
        srow.addWidget(self.add_substrate_btn)
        remove = QPushButton("Remove selected")
        remove.clicked.connect(self._remove_substrate)
        srow.addWidget(remove)
        srow.addStretch(1)
        bv.addLayout(srow)
        bv.addWidget(self._muted(
            "Starts from the model's main carbon and nitrogen sources. Add anything else "
            "worth tracking — iron, phosphate, a second C or N source. "
            "<b>Depletes</b>: consumed, and limits uptake as it runs out. "
            "<b>Held constant</b>: supplied in excess. Regulatory rules read these "
            "concentrations and the resulting fluxes automatically."))
        v.addWidget(box)

        run = QGroupBox("Culture")
        form = QFormLayout(run)
        self.tc_biomass = QDoubleSpinBox()
        self.tc_biomass.setRange(0.001, 10.0)
        self.tc_biomass.setDecimals(3)
        self.tc_biomass.setValue(0.05)
        self.tc_biomass.setMaximumWidth(110)
        form.addRow("Initial biomass (gDW L⁻¹):", self.tc_biomass)

        self.tc_duration = QDoubleSpinBox()
        self.tc_duration.setRange(1.0, 1000.0)
        self.tc_duration.setValue(96.0)
        self.tc_duration.setDecimals(0)
        self.tc_duration.setMaximumWidth(110)
        form.addRow("Duration (h):", self.tc_duration)

        self.tc_step = QDoubleSpinBox()
        self.tc_step.setRange(0.1, 24.0)
        self.tc_step.setValue(4.0)
        self.tc_step.setDecimals(1)
        self.tc_step.setMaximumWidth(110)
        self.tc_step.setToolTip("Smaller steps are more accurate and slower.")
        form.addRow("Step (h):", self.tc_step)

        self.tc_light = QDoubleSpinBox()
        self.tc_light.setRange(0.0, 3000.0)
        self.tc_light.setValue(300.0)
        self.tc_light.setDecimals(0)
        self.tc_light.setMaximumWidth(110)
        form.addRow("Light (µE):", self.tc_light)

        trow = QHBoxLayout()
        self.tc_targets_btn = QPushButton("Choose reactions…")
        self.tc_targets_btn.setToolTip(
            "Fluxes to record at every time point, besides growth and the tracked "
            "concentrations.")
        self.tc_targets_btn.clicked.connect(lambda: self._pick_targets("timecourse"))
        trow.addWidget(self.tc_targets_btn)
        self.tc_targets_label = QLabel("growth only")
        self.tc_targets_label.setStyleSheet(f"color:{style.TEXT_MUTED};")
        trow.addWidget(self.tc_targets_label, 1)
        form.addRow("Report:", trow)
        self._fix_height(self.tc_biomass, self.tc_duration, self.tc_step, self.tc_light,
                         self.tc_targets_btn, self.add_substrate_btn)
        v.addWidget(run)

        return self._page(
            page,
            "Explicit Euler — reliable for <i>when</i> a culture changes state, not for "
            "precise kinetics.")

    # -- stored-run bookkeeping ----------------------------------------------------
    def _unique_name(self, kind: str) -> str:
        stem = "ConditionScan" if kind == "scan" else "TimeCourse"
        self._counters[kind] = self._counters.get(kind, 0) + 1
        name = f"{stem}_{self._counters[kind]}"
        while name in self._runs:               # survives renames that collide
            self._counters[kind] += 1
            name = f"{stem}_{self._counters[kind]}"
        return name

    def add_run(self, frame, title: str, commentary: str = "", warning: str = "",
                kind: str = "") -> str:
        """Store a finished run under its own tab and select it."""
        name = self._unique_name("scan" if kind == "scan" else "timecourse")
        self._runs[name] = {"frame": frame, "title": title, "commentary": commentary,
                            "warning": warning, "kind": kind}
        index = self.run_tabs.addTab(name)
        self.run_tabs.setTabToolTip(index, title)
        self.run_tabs.setCurrentIndex(index)
        self._last_kind = kind
        self._sync_result_buttons()
        return name

    def _close_run(self, index: int) -> None:
        name = self.run_tabs.tabText(index)
        self._runs.pop(name, None)
        self.run_tabs.removeTab(index)
        self._sync_result_buttons()

    def _rename_run(self, index: int) -> None:
        if index < 0:
            return
        old = self.run_tabs.tabText(index)
        new, ok = QInputDialog.getText(self, "Rename run", "Name:", text=old)
        new = (new or "").strip()
        if not ok or not new or new == old:
            return
        if new in self._runs:
            QMessageBox.information(self, "Rename run",
                                    f"A run called “{new}” already exists.")
            return
        # Rebuilt rather than reassigned so the tab order and the record order stay the
        # same sequence; the plot dialog lists runs in tab order.
        self._runs = OrderedDict((new if key == old else key, value)
                                 for key, value in self._runs.items())
        self.run_tabs.setTabText(index, new)
        self._sync_result_buttons()

    def run_names(self) -> List[str]:
        return [self.run_tabs.tabText(i) for i in range(self.run_tabs.count())]

    def current_run(self) -> str:
        index = self.run_tabs.currentIndex()
        return self.run_tabs.tabText(index) if index >= 0 else ""

    def run_record(self, name: str) -> Optional[dict]:
        return self._runs.get(name)

    def runs_for_plot(self) -> List[tuple]:
        """(name, frame) for every stored run, in tab order."""
        return [(name, self._runs[name]["frame"]) for name in self.run_names()
                if name in self._runs]

    def _sync_result_buttons(self) -> None:
        has = bool(self._runs)
        for button in (getattr(self, "plot_btn", None), getattr(self, "table_btn", None)):
            if button is not None:
                button.setEnabled(has)
        hint = getattr(self, "run_hint", None)
        if hint is not None:
            hint.setText("" if has else
                         "Completed runs appear here. Double-click a tab to rename it.")

    # -- model awareness -----------------------------------------------------------
    def set_model(self, model) -> None:
        """Rebuild the condition list from this model's exchanges."""
        self._model = model
        self._sync_variables()
        self._clear_substrates()
        if model is not None:
            self._seed_substrates()

    def _exchanges(self) -> List[tuple]:
        """``(id, label, in_medium)`` for every exchange, medium components first."""
        if self._model is None:
            return []
        from ...core.network_graph import clean_label

        medium = set(getattr(self._model, "medium", {}) or {})
        rows = []
        for rxn in self._model.exchanges:
            label = clean_label(rxn.name or reg.exchange_metabolite_id(rxn.id))
            rows.append((rxn.id, label, rxn.id in medium))
        rows.sort(key=lambda r: (not r[2], r[1].lower()))
        return rows

    def _sync_variables(self) -> None:
        from ...core.network_graph import clean_label

        current = self.scan_variable.currentData()
        self.scan_variable.clear()
        for key, label, low, high in BASE_VARIABLES:
            self.scan_variable.addItem(label, {"kind": "sensor", "key": key,
                                               "low": low, "high": high})
        for rid, label, in_medium in self._exchanges():
            mark = "" if in_medium else "  (not in medium)"
            self.scan_variable.addItem(
                f"{label} — {clean_label(rid)}{mark}",
                {"kind": "exchange", "key": rid, "low": 0.0, "high": 50.0})
        index = self.scan_variable.findData(current)
        self.scan_variable.setCurrentIndex(max(0, index))
        self._sync_scan_range()

    def _sync_scan_range(self) -> None:
        data = self.scan_variable.currentData() or {}
        self.scan_from.setValue(float(data.get("low", 0.0)))
        self.scan_to.setValue(float(data.get("high", 100.0)))

    # -- substrate table -----------------------------------------------------------
    def _clear_substrates(self) -> None:
        self.substrates.setRowCount(0)

    def _seed_substrates(self) -> None:
        """Start with the model's main carbon and nitrogen source — nothing else.

        A batch culture is defined by what runs out, and in practice that is C or N. The
        two are chosen from the medium *this* model declares, taking the largest uptake
        allowance in each class, so the starting table describes the loaded model rather
        than a guess made here. Everything else — iron, phosphate, a second source — is
        one click away via <i>Add component…</i>, which is the right default: tracking a
        nutrient that never limits adds a column of flat numbers.
        """
        if self._model is None:
            return
        # Count carbon by parsing elements, not by looking for the letter "C": calcium
        # is "Ca", and a substring test picked EX_ca2_e as this model's carbon source.
        from ...core.physiology import _carbon_count

        medium = getattr(self._model, "medium", {}) or {}
        carbon, nitrogen = None, None
        carbon_cap = nitrogen_cap = -1.0

        for rxn in self._model.exchanges:
            allowance = float(medium.get(rxn.id, 0.0) or 0.0)
            if allowance <= 0:
                continue
            metabolite = next(iter(rxn.metabolites), None)
            if metabolite is None:
                continue
            sensor = reg.sensor_for_exchange(rxn.id)
            if sensor == "nitrogen_mM":
                if allowance > nitrogen_cap:
                    nitrogen, nitrogen_cap = rxn.id, allowance
            elif _carbon_count(metabolite) > 0 and allowance > carbon_cap:
                carbon, carbon_cap = rxn.id, allowance

        for exchange in (carbon, nitrogen):
            if exchange:
                self._add_substrate_row(exchange)

    def _add_substrate_row(self, exchange: str) -> None:
        if self._model is None or not self._model.reactions.has_id(exchange):
            return
        for row in range(self.substrates.rowCount()):
            cell = self.substrates.item(row, 1)
            if (cell.data(Qt.UserRole) or cell.text()) == exchange:
                return
        from ...core import screening as scr

        from ...core.network_graph import clean_label

        rxn = self._model.reactions.get_by_id(exchange)
        stem = reg.exchange_metabolite_id(exchange).rsplit("_", 1)[0]
        initial = DEFAULT_CONCENTRATION.get(stem, 10.0)
        # Shows the value the run will actually use, including the substitution made for
        # an unbounded medium entry — so the number in the table is the number applied.
        uptake = scr.default_uptake(self._model, exchange)

        row = self.substrates.rowCount()
        self.substrates.insertRow(row)
        # clean_label decodes the old BiGG `_LPAREN_e_RPAREN_` compartment encoding.
        # Without it, iJN678-derived entries read as "EX hco3 LPAREN e RPAREN".
        name = QTableWidgetItem(clean_label(rxn.name or stem))
        name.setFlags(name.flags() & ~Qt.ItemIsEditable)
        self.substrates.setItem(row, 0, name)
        rid = QTableWidgetItem(clean_label(exchange))
        rid.setFlags(rid.flags() & ~Qt.ItemIsEditable)
        rid.setData(Qt.UserRole, exchange)      # the real id, for the run
        self.substrates.setItem(row, 1, rid)
        self.substrates.setItem(row, 2, QTableWidgetItem(f"{initial:g}"))
        self.substrates.setItem(row, 3, QTableWidgetItem(f"{uptake:g}"))

        mode = QComboBox()
        mode.addItems([_MODE_DEPLETE, _MODE_BUFFERED])
        mode.setToolTip(
            "Depletes: the pool is consumed and limits uptake as it empties.\n"
            "Held constant: supplied in excess, so uptake keeps the medium's bound "
            "and the concentration does not fall.")
        self.substrates.setCellWidget(row, 4, mode)

    def _pick_substrates(self) -> None:
        if self._model is None:
            return
        from ..dialogs.reaction_browser_dialog import ReactionBrowserDialog

        existing = [self.substrates.item(r, 1).data(Qt.UserRole)
                    or self.substrates.item(r, 1).text()
                    for r in range(self.substrates.rowCount())]
        picked = ReactionBrowserDialog.pick(
            self, self._model, existing,
            title="Medium components to follow",
            prompt="Choose any exchange to follow through the run. Components already "
                   "in the medium are the usual choice, but any exchange can be "
                   "tracked.",
            restrict_to=[r.id for r in self._model.exchanges])
        if picked is None:
            return
        keep = set(picked)
        for row in reversed(range(self.substrates.rowCount())):
            cell = self.substrates.item(row, 1)
            if (cell.data(Qt.UserRole) or cell.text()) not in keep:
                self.substrates.removeRow(row)
        for rid in picked:
            self._add_substrate_row(rid)

    def _remove_substrate(self) -> None:
        for row in sorted({i.row() for i in self.substrates.selectedIndexes()},
                          reverse=True):
            self.substrates.removeRow(row)

    def substrate_specs(self) -> List[dict]:
        """What the time course should follow, read out of the table."""
        out = []
        for row in range(self.substrates.rowCount()):
            try:
                initial = float(self.substrates.item(row, 2).text())
            except (TypeError, ValueError):
                initial = 0.0
            try:
                uptake = float(self.substrates.item(row, 3).text())
            except (TypeError, ValueError):
                uptake = 10.0
            mode_widget = self.substrates.cellWidget(row, 4)
            cell = self.substrates.item(row, 1)
            out.append({
                # The shown id is cleaned for reading; the run needs the real one.
                "exchange": cell.data(Qt.UserRole) or cell.text(),
                "label": self.substrates.item(row, 0).text(),
                "initial_mM": initial,
                "max_uptake": uptake,
                # Left empty on purpose: which regulator reads a nutrient is a fact about
                # the organism, derived from the compound, not a user setting. Offering it
                # as a dropdown invited nonsense such as routing phosphate to the sulfur
                # sensor, which is what prompted removing the column.
                "sensor": "",
                "buffered": (mode_widget.currentText() == _MODE_BUFFERED
                             if mode_widget else False),
            })
        return out

    # -- target picking ------------------------------------------------------------
    def _pick_targets(self, which: str) -> None:
        if self._model is None:
            return
        from ..dialogs.reaction_browser_dialog import ReactionBrowserDialog

        current = self._scan_targets if which == "scan" else self._tc_targets
        picked = ReactionBrowserDialog.pick(
            self, self._model, current,
            title="Reactions to report",
            prompt="Growth is always reported. Choose any other fluxes to follow — "
                   "search by name, or filter to a subsystem and add it whole.")
        if picked is None:
            return
        if which == "scan":
            self._scan_targets = picked
            self._describe_targets(self.scan_targets_label, picked)
        else:
            self._tc_targets = picked
            self._describe_targets(self.tc_targets_label, picked)

    @staticmethod
    def _describe_targets(label: QLabel, targets: List[str]) -> None:
        if not targets:
            label.setText("growth only")
            label.setToolTip("")
            return
        shown = ", ".join(targets[:3])
        more = f" +{len(targets) - 3} more" if len(targets) > 3 else ""
        label.setText(f"growth + {len(targets)}: {shown}{more}")
        label.setToolTip("\n".join(targets))

    # -- helpers -------------------------------------------------------------------
    @staticmethod
    def _fix_height(*widgets) -> None:
        """Stop a control being compressed below its natural height.

        QFormLayout will happily shrink rows past their sizeHint when the container is
        short, which is what produced rows of half-visible text. Fixing the vertical
        policy makes the layout ask for more room (and the scroll area provide it)
        instead of clipping.
        """
        for widget in widgets:
            widget.setSizePolicy(widget.sizePolicy().horizontalPolicy(),
                                 QSizePolicy.Fixed)
            widget.setMinimumHeight(max(widget.sizeHint().height(), 24))

    def _muted(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setWordWrap(True)
        label.setTextFormat(Qt.RichText)
        label.setStyleSheet(f"color:{style.TEXT_MUTED};")
        return label

    def set_regulation_status(self, text: str, enabled: Optional[bool] = None) -> None:
        """Show the active rule set in one line; the full text lives in the tooltip."""
        plain = re.sub(r"<[^>]+>", "", text or "")
        first = plain.split(" — ")[0].split(".")[0].strip()
        if len(first) > 60:
            first = first[:57] + "…"
        self.rules_status.setText(f"<i>{first}</i>" if first else "")
        self.rules_status.setToolTip(plain)
        if enabled is not None:
            self.use_regulation.setChecked(bool(enabled))

    def regulation_enabled(self) -> bool:
        return bool(self.use_regulation.isChecked())

    # -- emitting ------------------------------------------------------------------
    def _sync_run_button(self) -> None:
        """One Run button, labelled with what it will actually run."""
        button = getattr(self, "run_btn", None)
        if button is None:
            return
        scan = self.controls.currentIndex() == 0
        button.setText("Run condition scan" if scan else "Run time course")

    def _emit_run(self) -> None:
        if self.controls.currentIndex() == 0:
            self._emit_scan()
        else:
            self._emit_timecourse()

    def _emit_plot(self) -> None:
        self.plot_requested.emit(self.run_names())

    def _emit_table(self) -> None:
        name = self.current_run()
        if name:
            self.table_requested.emit(name)

    def _emit_scan(self) -> None:
        data = self.scan_variable.currentData() or {}
        self.scan_requested.emit({
            "kind": data.get("kind", "sensor"),
            "variable": data.get("key", ""),
            "label": self.scan_variable.currentText(),
            "from": self.scan_from.value(),
            "to": self.scan_to.value(),
            "points": self.scan_points.value(),
            "logarithmic": bool(self.scan_log.currentData()),
            "targets": list(self._scan_targets),
            "regulation": self.regulation_enabled(),
        })

    def _emit_timecourse(self) -> None:
        self.timecourse_requested.emit({
            "substrates": self.substrate_specs(),
            "biomass": self.tc_biomass.value(),
            "duration_h": self.tc_duration.value(),
            "step_h": self.tc_step.value(),
            "light_uE": self.tc_light.value(),
            "targets": list(self._tc_targets),
            "regulation": self.regulation_enabled(),
        })

    # -- results -------------------------------------------------------------------
    def show_table(self, frame, title: str, commentary: str = "", warning: str = "",
                   kind: str = "") -> str:
        """Store a finished run. Kept under the old name so callers need not change."""
        return self.add_run(frame, title, commentary=commentary, warning=warning,
                            kind=kind)

    def set_busy(self, busy: bool) -> None:
        self.run_btn.setEnabled(not busy)
