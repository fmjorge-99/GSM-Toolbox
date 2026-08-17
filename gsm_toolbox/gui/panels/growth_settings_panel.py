"""Growth settings tab (Issue 5, redesigned in R4).

A dedicated, clearly-designed place to choose *what physiology the model is set
up for* before running any analysis — instead of silently inheriting whatever
medium/objective the SBML shipped with (which, for a cyanobacterium like iJN678,
is dark heterotrophy on glucose).

Layout (top → bottom):
  * live summary of the detected physiology + warnings,
  * biomass-reaction chooser (models often ship several),
  * the medium (editable uptake table) with a "Select Substrates…" shortcut that
    opens a dual-panel picker,
  * quick-mode presets at the BOTTOM — Autotrophic / Mixotrophic / Heterotrophic,
    each enabled only when that mode is possible for the loaded model.

It never edits the model directly; it emits signals the main window applies
through the undo-able project layer.
"""

from __future__ import annotations

import cobra
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ...core import media, physiology
from ..dialogs.substrate_picker_dialog import SubstratePickerDialog


class GrowthSettingsPanel(QWidget):
    # medium dict, biomass id (or None to keep current)
    apply_requested = Signal(dict, object)
    mode_requested = Signal(str)     # "autotrophic" | "mixotrophic" | "heterotrophic"

    def __init__(self):
        super().__init__()
        self._model = None
        self._spins = {}

        self.summary = QLabel("Load a model to configure its growth conditions.")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.summary.setStyleSheet(
            "background:#F5F6F8; border:1px solid #E0E3E8; border-radius:8px; padding:8px;")

        # ---- biomass objective ----
        bio_box = QGroupBox("Biomass / growth objective")
        bio_layout = QHBoxLayout(bio_box)
        bio_layout.addWidget(QLabel("Objective:"))
        self.biomass_combo = QComboBox()
        self.biomass_combo.setMinimumWidth(280)
        bio_layout.addWidget(self.biomass_combo, 1)
        self.set_biomass_btn = QPushButton("Set objective")
        self.set_biomass_btn.clicked.connect(self._emit_apply)
        bio_layout.addWidget(self.set_biomass_btn)

        # ---- medium table + substrate picker ----
        med_box = QGroupBox("Medium — maximum uptake (mmol/gDW/h; 0 = not supplied)")
        med_layout = QVBoxLayout(med_box)
        picker_row = QHBoxLayout()
        self.pick_btn = QPushButton("Select Substrates…")
        self.pick_btn.setToolTip("Quickly choose which nutrients the cell may take up "
                                 "(a shortcut; fine-tune rates in the table below).")
        self.pick_btn.clicked.connect(self._open_substrate_picker)
        picker_row.addWidget(self.pick_btn)
        picker_row.addStretch(1)
        med_layout.addLayout(picker_row)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Exchange", "Name", "Max uptake"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        med_layout.addWidget(self.table, 1)
        apply_row = QHBoxLayout()
        apply_row.addStretch(1)
        self.apply_btn = QPushButton("Apply medium")
        self.apply_btn.setToolTip("Apply the uptake bounds above as the growth medium.")
        self.apply_btn.clicked.connect(self._emit_apply)
        apply_row.addWidget(self.apply_btn)
        med_layout.addLayout(apply_row)

        # ---- quick-mode presets (bottom) ----
        preset_box = QGroupBox("Quick presets — growth mode")
        preset_layout = QHBoxLayout(preset_box)
        self.auto_btn = QPushButton("Autotrophic")
        self.auto_btn.setToolTip("Light + inorganic carbon (CO₂/HCO₃⁻), autotrophic biomass.")
        self.auto_btn.clicked.connect(lambda: self.mode_requested.emit("autotrophic"))
        self.mixo_btn = QPushButton("Mixotrophic")
        self.mixo_btn.setToolTip("Light AND an organic carbon source together.")
        self.mixo_btn.clicked.connect(lambda: self.mode_requested.emit("mixotrophic"))
        self.hetero_btn = QPushButton("Heterotrophic")
        self.hetero_btn.setToolTip("Organic carbon only, no light (dark growth).")
        self.hetero_btn.clicked.connect(lambda: self.mode_requested.emit("heterotrophic"))
        for b in (self.auto_btn, self.mixo_btn, self.hetero_btn):
            b.setMinimumHeight(36)
            preset_layout.addWidget(b)

        inner = QWidget()
        v = QVBoxLayout(inner)
        v.addWidget(self.summary)
        v.addWidget(bio_box)
        v.addWidget(med_box, 1)
        v.addWidget(preset_box)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(inner)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)

        self._set_enabled(False)

    def _set_enabled(self, on: bool) -> None:
        for w in (self.set_biomass_btn, self.apply_btn, self.biomass_combo, self.table,
                  self.pick_btn, self.auto_btn, self.mixo_btn, self.hetero_btn):
            w.setEnabled(on)

    def set_model(self, model: cobra.Model) -> None:
        self._model = model
        if model is None:
            self._set_enabled(False)
            self.summary.setText("Load a model to configure its growth conditions.")
            return
        self._set_enabled(True)
        self._refresh_summary()
        self._populate_biomass()
        self._populate_medium()
        self._update_mode_availability()

    def _refresh_summary(self) -> None:
        s = physiology.summarize(self._model)
        text = (f"<b>Detected physiology</b><br>Objective: <b>{s.objective}</b><br>"
                f"Energy source: <b>{s.energy_source}</b> · Carbon: <b>{s.carbon_source}</b>")
        if s.warnings:
            text += "<br>" + "<br>".join(f"⚠ {w}" for w in s.warnings)
        self.summary.setText(text)

    def _update_mode_availability(self) -> None:
        modes = physiology.available_growth_modes(self._model)
        self.auto_btn.setEnabled(modes["autotrophic"])
        self.mixo_btn.setEnabled(modes["mixotrophic"])
        self.hetero_btn.setEnabled(modes["heterotrophic"])
        if not modes["autotrophic"]:
            self.auto_btn.setToolTip("No photon exchange — this model is not phototroph-capable.")
        if not modes["mixotrophic"]:
            self.mixo_btn.setToolTip("Mixotrophy needs both light and an organic carbon source.")
        if not modes["heterotrophic"]:
            self.hetero_btn.setToolTip("No organic carbon exchange found in this model.")

    def _populate_biomass(self) -> None:
        self.biomass_combo.blockSignals(True)
        self.biomass_combo.clear()
        self.biomass_combo.addItem("(keep current objective)", userData=None)
        for b in physiology.list_biomass_reactions(self._model):
            tag = f" [{b['kind']}]" if b["kind"] else ""
            cur = " (current)" if b["objective"] else ""
            label = f"{b['name'] or b['id']} ({b['id']}){tag}{cur}"
            self.biomass_combo.addItem(label, userData=b["id"])
        self.biomass_combo.blockSignals(False)

    def _populate_medium(self) -> None:
        rows = media.list_exchanges(self._model)
        rows.sort(key=lambda r: (r["uptake"] <= 0, r["id"]))   # active uptakes first
        self.table.setRowCount(len(rows))
        self._spins = {}
        for r, row in enumerate(rows):
            self.table.setItem(r, 0, _ro_item(row["id"]))
            self.table.setItem(r, 1, _ro_item(row["name"]))
            spin = QDoubleSpinBox()
            spin.setRange(0.0, 1e6)
            spin.setDecimals(3)
            spin.setValue(float(row["uptake"]))
            self._spins[row["id"]] = spin
            self.table.setCellWidget(r, 2, spin)

    def _open_substrate_picker(self) -> None:
        if self._model is None:
            return
        dlg = SubstratePickerDialog(self._model, self)
        if dlg.exec() != SubstratePickerDialog.Accepted:
            return
        chosen = set(dlg.selected_ids())
        default_uptake = 10.0
        # Turn selected exchanges on (default rate if currently off), others off.
        for rid, spin in self._spins.items():
            if rid in chosen:
                if spin.value() <= 0:
                    spin.setValue(default_uptake)
            else:
                spin.setValue(0.0)
        self._emit_apply()

    def _medium(self) -> dict:
        return {rid: s.value() for rid, s in self._spins.items() if s.value() > 0}

    def _emit_apply(self) -> None:
        self.apply_requested.emit(self._medium(), self.biomass_combo.currentData())


def _ro_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item
