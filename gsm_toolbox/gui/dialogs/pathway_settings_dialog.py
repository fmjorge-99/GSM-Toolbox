"""Settings dialog shown when predicting a heterologous pathway.

Lets the user control the network-expansion search: which metabolites to start
from, preferred enzyme classes, how long a pathway may be, and how many
alternative routes to return.
"""

from __future__ import annotations

from typing import List, Tuple

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


class PathwaySettingsDialog(QDialog):
    """Collect pathway-prediction options before running the search."""

    def __init__(self, parent, target_display: str, metabolites: List[Tuple[str, str]],
                 *, min_flux: float = 0.1):
        """``metabolites`` is a list of ``(display, id)`` pairs from the model."""
        super().__init__(parent)
        self.setWindowTitle("Pathway design settings")
        self.setMinimumWidth(520)
        self._metabolites = metabolites

        layout = QVBoxLayout(self)

        intro = QLabel(
            f"Predict a heterologous pathway to <b>{target_display}</b>. By default the search "
            "starts from everything your model already makes and finds the shortest route. Adjust "
            "the options below for finer control.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # --- search options FIRST (#B10) -----------------------------------
        form_box = QGroupBox("1 · Search options")
        form = QFormLayout(form_box)

        self.algorithm_combo = QComboBox()
        self.algorithm_combo.addItem("Retrosynthesis (target → native precursors)", "retro")
        self.algorithm_combo.addItem("Network expansion (host → target)", "expansion")
        self.algorithm_combo.setToolTip(
            "Retrosynthesis works backward from the target to what the host already "
            "makes, giving the minimal set of heterologous steps (recommended). "
            "Network expansion grows forward from the host's metabolites.")
        form.addRow("Search algorithm:", self.algorithm_combo)

        self.n_alts = QSpinBox()
        self.n_alts.setRange(1, 20)
        self.n_alts.setValue(1)
        self.n_alts.setToolTip(
            "Return this many distinct alternative pathways. After finding the best route, its "
            "reactions are excluded and the search repeats to yield genuinely different options.")
        self.n_alts.valueChanged.connect(self._update_alt_options)
        form.addRow("Number of pathways to return:", self.n_alts)

        self.priority_combo = QComboBox()
        self.priority_combo.addItem("Highest predicted flux (yield)", "yield")
        self.priority_combo.addItem("Fewest reactions (shortest)", "shortest")
        self.priority_combo.setToolTip("How to rank alternative routes against each other "
                                       "(only when returning more than one pathway).")
        self._priority_label = QLabel("Prioritise:")
        form.addRow(self._priority_label, self.priority_combo)

        self.max_steps = QSpinBox()
        self.max_steps.setRange(1, 60)
        self.max_steps.setValue(25)
        self.max_steps.setToolTip("Maximum number of heterologous reactions in a pathway. Larger "
                                  "values explore longer routes but take more time.")
        form.addRow("Maximum pathway length:", self.max_steps)

        self.include_boundary_check = QCheckBox(
            "Allow database exchange/demand reactions in routes")
        self.include_boundary_check.setToolTip(
            "Off (recommended): routes use only real metabolic reactions. On: the database's own "
            "exchange/demand reactions may be used, which can let a metabolite appear from an "
            "external source.")
        form.addRow("Boundary reactions:", self.include_boundary_check)

        self.ec_edit = QLineEdit()
        self.ec_edit.setPlaceholderText("e.g. 1.1.1, 2.3.1.- (comma-separated)")
        self.ec_edit.setToolTip(
            "Preferred EC classes. Reactions whose enzyme classification matches are tried first, "
            "biasing the route toward chemistries you favour. Leave blank for no preference.")
        form.addRow("Preferred EC classes:", self.ec_edit)

        self.min_flux = QDoubleSpinBox()
        self.min_flux.setRange(0.0, 1000.0)
        self.min_flux.setDecimals(3)
        self.min_flux.setValue(min_flux)
        self.min_flux.setToolTip("Target production flux the pathway should enable (informational "
                                 "— the predicted flux is reported after the search).")
        form.addRow("Minimum target flux:", self.min_flux)
        layout.addWidget(form_box)

        # --- starting metabolites SECOND -----------------------------------
        start_box = QGroupBox("2 · Starting metabolites (optional)")
        start_v = QVBoxLayout(start_box)
        self.use_all_check = QCheckBox("Use all metabolites/reactions in the model (default)")
        self.use_all_check.setChecked(True)
        self.use_all_check.setToolTip(
            "Start the search from every metabolite the model can make. Uncheck to restrict the "
            "pathway to originate from a specific set of precursors you choose below.")
        self.use_all_check.toggled.connect(self._toggle_start)
        start_v.addWidget(self.use_all_check)

        # A model contains many compounds that exist in the stoichiometry but carry no
        # flux while the cell grows. Letting a route "start" from one of those produces a
        # design that assumes a precursor the host does not actually make — lactaldehyde
        # in Synechocystis is the classic case, and without this option the search simply
        # reports that 1,2-propanediol needs no heterologous reactions at all.
        self.flux_only_check = QCheckBox(
            "Only start from metabolites that carry flux in growth conditions")
        self.flux_only_check.setChecked(False)
        self.flux_only_check.setToolTip(
            "Restrict starting points to compounds the host actually produces at the "
            "current growth solution.\n\n"
            "Many metabolites are present in a model but idle — nothing makes them under "
            "these conditions. A route that begins at one of those is not buildable "
            "without extra genes.\n\n"
            "With this on, the search must bridge the gap itself, adding the "
            "heterologous steps that connect an available precursor to the target.")
        start_v.addWidget(self.flux_only_check)
        flux_note = QLabel(
            "Costs one extra FBA solve. Leave off to reproduce earlier results, or when "
            "your model has no meaningful growth objective.")
        flux_note.setWordWrap(True)
        flux_note.setStyleSheet("color:#5f6368;")
        start_v.addWidget(flux_note)

        picker_row = QHBoxLayout()
        self.met_combo = QComboBox()
        self.met_combo.setEditable(True)
        self.met_combo.setInsertPolicy(QComboBox.NoInsert)
        for display, mid in metabolites:
            self.met_combo.addItem(display, mid)
        self.add_met_btn = QPushButton("Add →")
        self.add_met_btn.clicked.connect(self._add_metabolite)
        picker_row.addWidget(self.met_combo, 1)
        picker_row.addWidget(self.add_met_btn)
        start_v.addLayout(picker_row)

        self.start_list = QListWidget()
        self.start_list.setMaximumHeight(120)
        start_v.addWidget(self.start_list)
        remove_btn = QPushButton("Remove selected")
        remove_btn.clicked.connect(self._remove_selected)
        start_v.addWidget(remove_btn)
        self._start_widgets = [self.met_combo, self.add_met_btn, self.start_list, remove_btn]
        layout.addWidget(start_box)
        self._update_alt_options()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        buttons.button(QDialogButtonBox.Ok).setText("Predict")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._toggle_start(True)

    def _toggle_start(self, use_all: bool) -> None:
        for w in self._start_widgets:
            w.setEnabled(not use_all)

    def _update_alt_options(self) -> None:
        # "Prioritise" only matters when returning more than one pathway (#B10).
        multi = self.n_alts.value() > 1
        self.priority_combo.setEnabled(multi)
        self._priority_label.setEnabled(multi)

    def _add_metabolite(self) -> None:
        mid = self.met_combo.currentData()
        display = self.met_combo.currentText()
        if mid is None:
            # user typed a raw id
            mid = display.strip()
            if not mid:
                return
        # avoid duplicates
        for i in range(self.start_list.count()):
            if self.start_list.item(i).data(256) == mid:
                return
        item_text = display if display else str(mid)
        self.start_list.addItem(item_text)
        self.start_list.item(self.start_list.count() - 1).setData(256, mid)

    def _remove_selected(self) -> None:
        for item in self.start_list.selectedItems():
            self.start_list.takeItem(self.start_list.row(item))

    def values(self) -> dict:
        start_mets: List[str] = []
        if not self.use_all_check.isChecked():
            for i in range(self.start_list.count()):
                start_mets.append(str(self.start_list.item(i).data(256)))
        ec = [e.strip() for e in self.ec_edit.text().split(",") if e.strip()]
        return {
            "use_all": self.use_all_check.isChecked(),
            "start_metabolites": start_mets,
            "flux_carrying_starts_only": self.flux_only_check.isChecked(),
            "preferred_ec": ec,
            "max_steps": self.max_steps.value(),
            "n_alternatives": self.n_alts.value(),
            "min_flux": self.min_flux.value(),
            "algorithm": self.algorithm_combo.currentData(),
            "priority": self.priority_combo.currentData(),
            "include_boundary": self.include_boundary_check.isChecked(),
        }
