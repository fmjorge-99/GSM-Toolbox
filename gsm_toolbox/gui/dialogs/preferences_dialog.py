"""Settings ▸ Preferences — every configurable value in the toolbox, in one place.

The previous version was a single scrolling column of three option groups. As the app
grew that stopped scaling: options had no obvious home, and the dialog ran off the bottom
of the screen. This is a **categorised** dialog — a list of categories on the left, one
page each on the right — so a setting can be found by guessing its category rather than by
scrolling past everything else.

Categories are chosen so that every configurable value has exactly one plausible home:

* **General** — where results go, what happens on startup and exit
* **Appearance** — font, panel layout, table density
* **Pathway Finder** — everything the route search and its add-ons use
* **Analysis & Solver** — how the solver is run
* **Regulation & Dynamics** — the regulatory layer and time-course defaults
* **Data & Storage** — caches, offline mode, download consent
* **Shortcuts** — what appears on the quick-access bar

Long explanations stay in collapsible sections *within* a page, so a page can be as
informative as it needs to be without the dialog outgrowing the display.
"""
from __future__ import annotations

from typing import Dict, List

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMessageBox, QPushButton, QScrollArea, QSpinBox, QSplitter,
    QStackedWidget, QVBoxLayout, QWidget)

from ...core import preferences as prefs
from ...core import thermodynamics as thermo
from .. import style
from ..widgets.collapsible import CollapsibleSection

#: Actions offered for the quick-access bar. The id is what gets stored; the label is
#: what the user sees. Chosen as the things a user reaches for repeatedly rather than
#: everything that exists — a configurable bar of fifty buttons is not configurable.
TOOLBAR_CANDIDATES: List[tuple] = [
    ("open_model", "Open Model"),
    ("save_project", "Save Project"),
    ("add_reaction", "Add Reaction"),
    ("undo", "Undo"),
    ("redo", "Redo"),
    ("growth_settings", "Growth Settings"),
    ("network_visualization", "Network Visualization"),
    ("run_fba", "Run FBA"),
    ("run_pfba", "Run pFBA"),
    ("pathway_design", "Pathway Design"),
    ("dynamic_analysis", "Dynamic Analysis"),
    ("edit_medium", "Edit Growth Medium"),
    ("manage_databases", "Manage Databases"),
    ("preferences", "Preferences"),
]


def _muted(text: str) -> QLabel:
    label = QLabel(text)
    label.setWordWrap(True)
    label.setTextFormat(Qt.RichText)
    label.setStyleSheet(f"color:{style.TEXT_MUTED};")
    return label


class PreferencesDialog(QDialog):
    """Categorised settings. Nothing is written until OK is pressed."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.resize(940, 620)
        from ..widgets.dialog_util import clamp_to_screen

        self._widgets: Dict[str, object] = {}

        outer = QVBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)

        self.categories = QListWidget()
        self.categories.setMaximumWidth(220)
        self.categories.setSelectionMode(QAbstractItemView.SingleSelection)
        splitter.addWidget(self.categories)

        self.pages = QStackedWidget()
        splitter.addWidget(self.pages)
        splitter.setStretchFactor(1, 1)
        outer.addWidget(splitter, 1)

        for title, builder in (
                ("General", self._page_general),
                ("Appearance", self._page_appearance),
                ("Pathway Finder", self._page_pathway),
                ("Analysis & Solver", self._page_analysis),
                ("Regulation & Dynamics", self._page_regulation),
                ("Data & Storage", self._page_data),
                ("Shortcuts", self._page_shortcuts)):
            self.categories.addItem(QListWidgetItem(title))
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QScrollArea.NoFrame)
            page = QWidget()
            builder(QVBoxLayout(page))
            scroll.setWidget(page)
            self.pages.addWidget(scroll)

        self.categories.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.categories.setCurrentRow(0)

        # Outside the splitter, so no page length can push them off-screen.
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel
                                   | QDialogButtonBox.RestoreDefaults)
        buttons.accepted.connect(self._save)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self._restore)
        outer.addWidget(buttons)

        self._load()
        clamp_to_screen(self, frac=0.9)

    # -- small builders ------------------------------------------------------------
    def _check(self, key: str, label: str, tip: str = "") -> QCheckBox:
        box = QCheckBox(label)
        if tip:
            box.setToolTip(tip)
        self._widgets[key] = box
        return box

    def _spin(self, key: str, low: int, high: int, tip: str = "") -> QSpinBox:
        box = QSpinBox()
        box.setRange(low, high)
        if tip:
            box.setToolTip(tip)
        self._widgets[key] = box
        return box

    def _dspin(self, key: str, low: float, high: float, step: float, decimals: int,
               tip: str = "") -> QDoubleSpinBox:
        box = QDoubleSpinBox()
        box.setRange(low, high)
        box.setSingleStep(step)
        box.setDecimals(decimals)
        if tip:
            box.setToolTip(tip)
        self._widgets[key] = box
        return box

    def _combo(self, key: str, options: List[tuple], tip: str = "") -> QComboBox:
        box = QComboBox()
        for value, label in options:
            box.addItem(label, value)
        if tip:
            box.setToolTip(tip)
        self._widgets[key] = box
        return box

    # -- pages ---------------------------------------------------------------------
    def _page_general(self, v: QVBoxLayout) -> None:
        form = QFormLayout()
        row = QHBoxLayout()
        self.results_dir = QLineEdit()
        self.results_dir.setPlaceholderText("(not set — results are not saved automatically)")
        self._widgets[prefs.RESULTS_DIR] = self.results_dir
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._pick_results_dir)
        row.addWidget(self.results_dir, 1)
        row.addWidget(browse)
        form.addRow("Results output folder:", row)
        form.addRow("", self._check(
            prefs.AUTOSAVE_RESULTS, "Save every analysis table as CSV automatically",
            "Each analysis writes a timestamped CSV into the folder above."))
        form.addRow("", self._check(
            prefs.CONFIRM_ON_EXIT, "Ask before closing with unsaved changes"))
        form.addRow("", self._check(
            prefs.RESTORE_LAST_PROJECT, "Reopen the last project on startup"))
        v.addLayout(form)
        v.addWidget(_muted(
            "The results folder is where analysis tables are written. Leaving it unset "
            "does not lose anything — results stay in the app and can be saved by hand."))
        v.addStretch(1)

    def _page_appearance(self, v: QVBoxLayout) -> None:
        form = QFormLayout()
        form.addRow("Font size:", self._spin(
            prefs.FONT_SIZE, 7, 18, "Base point size for the whole interface."))
        form.addRow("Information panel:", self._combo(
            prefs.INFO_PANEL_POSITION,
            [("right", "Right side"), ("bottom", "Bottom (under centre)")],
            "Where the Information dock is placed. Moved here from the View menu."))
        form.addRow("Table density:", self._combo(
            prefs.TABLE_DENSITY,
            [("comfortable", "Comfortable"), ("compact", "Compact")],
            "Compact fits more rows on screen."))
        v.addLayout(form)

        panels = CollapsibleSection("Panels shown by default",
                                    summary="Explorer · Categories · Information")
        panels.add_widget(self._check(prefs.SHOW_EXPLORER, "Explorer"), always_visible=True)
        panels.add_widget(self._check(prefs.SHOW_CATEGORIES, "Categories"),
                          always_visible=True)
        panels.add_widget(self._check(prefs.SHOW_INFO, "Information"), always_visible=True)
        panels.add_widget(_muted(
            "These set what is visible when the app starts. Panels can still be shown or "
            "hidden at any time from the View menu."))
        v.addWidget(panels)
        v.addStretch(1)

    def _page_pathway(self, v: QVBoxLayout) -> None:
        search = CollapsibleSection("Search defaults",
                                    summary="algorithm, depth, alternatives")
        form = QFormLayout()
        form.addRow("Algorithm:", self._combo(
            prefs.SEARCH_ALGORITHM,
            [("retro", "Retrosynthetic (default)"), ("expansion", "Forward expansion")]))
        form.addRow("Maximum steps:", self._spin(prefs.SEARCH_MAX_STEPS, 1, 100))
        form.addRow("Alternatives to return:", self._spin(prefs.SEARCH_ALTERNATIVES, 1, 20))
        search.add_layout(form, always_visible=True)
        search.add_widget(self._check(
            prefs.FLUX_CARRYING_STARTS,
            "Only start from metabolites that carry flux in growth conditions",
            "Stops a route being anchored on a compound the host does not actually make."))
        search.add_widget(_muted(
            "A model can contain compounds it never produces. Restricting starting points "
            "to flux-carrying metabolites forces the search to include the heterologous "
            "steps that supply them — but changes results, so it is off by default."))
        v.addWidget(search)

        rr = CollapsibleSection("RetroRules", summary="rule-based retrosynthesis")
        rrform = QFormLayout()
        rrform.addRow("Search seed:", self._spin(
            prefs.RETRORULES_SEED, 0, 2_000_000_000,
            "Fixes the order rules are explored so a query returns the same routes "
            "every time."))
        rrform.addRow("Rule diameter:", self._spin(
            prefs.RR_DIAMETER, 2, 16,
            "Larger = more specific rules, fewer but more reliable matches."))
        rr.add_layout(rrform, always_visible=True)
        rr.add_widget(self._check(
            prefs.RR_PLAUSIBILITY, "Hide chemically implausible routes"))
        rr.add_widget(_muted(
            "The rule set is large and the search can only sample it. A fixed seed makes "
            "that sample reproducible; it does not make it exhaustive."))
        v.addWidget(rr)

        mdf = CollapsibleSection("Thermodynamics (MDF)",
                                 summary="off by default · one-off ~1.34 GB download")
        self.mdf_check = self._check(prefs.MDF_ENABLED, "Enable the MDF suite")
        mdf.add_widget(self.mdf_check, always_visible=True)
        self.mdf_status = QLabel("")
        self.mdf_status.setWordWrap(True)
        mdf.add_widget(self.mdf_status, always_visible=True)
        mdf.add_widget(_muted(
            "Estimates whether a route can be driven forward. Its answers depend on "
            "assumed pH, ionic strength and metabolite concentrations, so treat it as an "
            "advanced check rather than a headline number."))
        v.addWidget(mdf)

        sz = CollapsibleSection("Enzyme search (Selenzyme-style)",
                                summary="finds enzymes for reactions with no EC number")
        self.sz_check = self._check(prefs.SELENZYME_ENABLED,
                                    "Enable reaction-similarity enzyme search")
        sz.add_widget(self.sz_check, always_visible=True)
        self.sz_status = QLabel("")
        self.sz_status.setWordWrap(True)
        sz.add_widget(self.sz_status, always_visible=True)
        sz.add_widget(self._check(
            prefs.STRUCTURE_FETCH, "Fetch missing structures from the web when needed"))
        v.addWidget(sz)
        v.addStretch(1)

    def _page_analysis(self, v: QVBoxLayout) -> None:
        form = QFormLayout()
        form.addRow("Solver:", self._combo(
            prefs.SOLVER, [("", "Automatic"), ("glpk", "GLPK"), ("cplex", "CPLEX"),
                           ("gurobi", "Gurobi"), ("scipy", "SciPy")],
            "Leave on Automatic unless a specific solver is required."))
        form.addRow("FVA fraction of optimum:", self._dspin(
            prefs.FVA_FRACTION, 0.0, 1.0, 0.05, 2,
            "How much of the optimum must be retained during flux variability analysis."))
        form.addRow("Worker processes:", self._spin(
            prefs.WORKER_PROCESSES, 1, 32,
            "Parallel processes for FVA and similar. Values above 1 can be unreliable "
            "on Windows, which is why the default is 1."))
        form.addRow("Growth floor for production:", self._dspin(
            prefs.PRODUCTION_GROWTH_FLOOR, 0.0, 1.0, 0.05, 2,
            "Fraction of maximum growth a strain must retain when a product maximum is "
            "computed. Zero would report a number no living strain can hold."))
        v.addLayout(form)
        v.addWidget(_muted(
            "Maximum product fluxes are reported with growth held at this fraction of the "
            "condition's own maximum. Raising it gives more conservative, more realistic "
            "titres."))
        v.addStretch(1)

    def _page_regulation(self, v: QVBoxLayout) -> None:
        section = CollapsibleSection(
            "Regulatory layer",
            summary="off by default · no rule set until you load one")
        section.add_widget(self._check(
            prefs.REGULATION_ENABLED, "Apply regulatory rules during simulation",
            "Rules gate or scale reactions, enzyme costs and the protein budget as "
            "conditions change."), always_visible=True)
        row = QHBoxLayout()
        self.ruleset_path = QLineEdit()
        self.ruleset_path.setPlaceholderText("(no rule set active)")
        self._widgets[prefs.REGULATION_RULESET] = self.ruleset_path
        pick = QPushButton("Browse…")
        pick.clicked.connect(self._pick_ruleset)
        row.addWidget(self.ruleset_path, 1)
        row.addWidget(pick)
        section.add_layout(row, always_visible=True)
        section.add_widget(_muted(
            "Rule sets are files you load — nothing is applied until you choose one, "
            "because regulation describes a particular organism. Manage them from "
            "<i>Tools ▸ Regulation</i>. A result that depends on an <i>assumed</i> "
            "threshold is flagged wherever it is reported."))
        v.addWidget(section)

        form = QFormLayout()
        form.addRow("Rule activation threshold:", self._dspin(
            prefs.REGULATION_ACTIVATION, 0.0, 1.0, 0.01, 2,
            "Below this a rule is not counted as firing. A de-repression response never "
            "reaches exactly zero, so a threshold of 0 would mark every such rule as "
            "permanently active and hide real transitions."))
        form.addRow("Time-course step (h):", self._dspin(
            prefs.DFBA_STEP_H, 0.1, 24.0, 0.5, 1,
            "Integration step for dynamic simulations. Smaller is more accurate and "
            "slower."))
        form.addRow("Time-course duration (h):", self._dspin(
            prefs.DFBA_DURATION_H, 1.0, 1000.0, 12.0, 0))
        v.addLayout(form)
        v.addWidget(_muted(
            "Dynamic simulation uses explicit Euler integration — suitable for phase "
            "behaviour (when a culture switches state), not for precise kinetics."))
        v.addStretch(1)

    def _page_data(self, v: QVBoxLayout) -> None:
        v.addWidget(self._check(
            prefs.OFFLINE_MODE, "Offline mode — never contact the network",
            "Blocks structure fetching and reference-data downloads. Anything already "
            "cached still works."))
        v.addWidget(self._check(
            prefs.ALLOW_DOWNLOADS, "Allow one-off reference-data downloads",
            "Some optional features need a large dataset the first time they are used."))
        v.addWidget(_muted(
            "Cached databases and molecule images are managed from "
            "<b>Settings ▸ Manage Data…</b>, where they can be inspected and cleared."))
        v.addStretch(1)

    def _page_shortcuts(self, v: QVBoxLayout) -> None:
        v.addWidget(QLabel("<b>Quick-access bar</b>"))
        v.addWidget(_muted(
            "Tick the actions to show on the bar below the menus. The order follows the "
            "list."))
        self.shortcut_boxes: Dict[str, QCheckBox] = {}
        for action_id, label in TOOLBAR_CANDIDATES:
            box = QCheckBox(label)
            self.shortcut_boxes[action_id] = box
            v.addWidget(box)
        reset = QPushButton("Reset to the default bar")
        reset.clicked.connect(self._reset_toolbar)
        v.addWidget(reset)
        v.addStretch(1)

    # -- helpers -------------------------------------------------------------------
    def _pick_results_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Results output folder",
                                                self.results_dir.text() or "")
        if path:
            self.results_dir.setText(path)

    def _pick_ruleset(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Regulatory rule set", "",
                                              "Rule sets (*.json);;All files (*)")
        if path:
            self.ruleset_path.setText(path)

    def _reset_toolbar(self) -> None:
        for action_id, box in self.shortcut_boxes.items():
            box.setChecked(action_id in prefs.TOOLBAR_DEFAULT)

    def _refresh_status(self) -> None:
        from ...core import selenzyme as sz
        if thermo.cache_present():
            self.mdf_status.setText(
                "<span style='color:#188038'>✓ Thermodynamics data is installed.</span>")
        else:
            self.mdf_status.setText(
                "<span style='color:#E8710A'>Data not downloaded yet (~1.34 GB). It is "
                "fetched when the suite is first used.</span>")
        try:
            if not sz.rdkit_available():
                self.sz_status.setText(
                    "<span style='color:#C5221F'>RDKit is unavailable in this build, so "
                    "reaction-similarity search cannot be enabled.</span>")
                self.sz_check.setEnabled(False)
            elif sz.is_installed():
                self.sz_status.setText(
                    "<span style='color:#188038'>✓ Enzyme index is built.</span>")
            else:
                self.sz_status.setText(
                    "<span style='color:#E8710A'>Enzyme index not built yet; it is built "
                    "on first use.</span>")
        except Exception:  # noqa: BLE001 - status is informational only
            self.sz_status.setText("")

    # -- load / save ---------------------------------------------------------------
    def _load(self) -> None:
        for key, widget in self._widgets.items():
            value = prefs.get(key)
            if isinstance(widget, QCheckBox):
                widget.setChecked(bool(value))
            elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                widget.setValue(type(widget.value())(value))
            elif isinstance(widget, QComboBox):
                index = widget.findData(value)
                widget.setCurrentIndex(index if index >= 0 else 0)
            elif isinstance(widget, QLineEdit):
                widget.setText(str(value or ""))
        active = prefs.get(prefs.TOOLBAR_ACTIONS) or list(prefs.TOOLBAR_DEFAULT)
        for action_id, box in self.shortcut_boxes.items():
            box.setChecked(action_id in active)
        self._refresh_status()

    def _restore(self) -> None:
        if QMessageBox.question(
                self, "Restore defaults",
                "Reset every preference to its default value?") != QMessageBox.Yes:
            return
        for key, widget in self._widgets.items():
            value = prefs._DEFAULTS.get(key)
            if isinstance(widget, QCheckBox):
                widget.setChecked(bool(value))
            elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                widget.setValue(type(widget.value())(value))
            elif isinstance(widget, QComboBox):
                index = widget.findData(value)
                widget.setCurrentIndex(index if index >= 0 else 0)
            elif isinstance(widget, QLineEdit):
                widget.setText(str(value or ""))
        self._reset_toolbar()

    def _save(self) -> None:
        from ...core import selenzyme as sz
        for key, widget in self._widgets.items():
            if isinstance(widget, QCheckBox):
                prefs.set(key, bool(widget.isChecked()))
            elif isinstance(widget, QSpinBox):
                prefs.set(key, int(widget.value()))
            elif isinstance(widget, QDoubleSpinBox):
                prefs.set(key, float(widget.value()))
            elif isinstance(widget, QComboBox):
                prefs.set(key, widget.currentData())
            elif isinstance(widget, QLineEdit):
                prefs.set(key, widget.text().strip())

        # An optional suite stays off until its data exists, whatever the box says —
        # enabling it without the data would surface controls that cannot work.
        if prefs.get(prefs.MDF_ENABLED) and not thermo.cache_present():
            prefs.set(prefs.MDF_ENABLED, False)
            QMessageBox.information(
                self, "Thermodynamics data missing",
                "The MDF suite stays disabled until its dataset has been downloaded.")
        try:
            if prefs.get(prefs.SELENZYME_ENABLED) and not sz.is_installed():
                prefs.set(prefs.SELENZYME_ENABLED, False)
                QMessageBox.information(
                    self, "Enzyme index missing",
                    "Reaction-similarity search stays disabled until its index is built.")
        except Exception:  # noqa: BLE001
            pass

        chosen = [a for a, box in self.shortcut_boxes.items() if box.isChecked()]
        prefs.set(prefs.TOOLBAR_ACTIONS, chosen or list(prefs.TOOLBAR_DEFAULT))
        self.accept()
