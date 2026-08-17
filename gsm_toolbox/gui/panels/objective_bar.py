"""A user-friendly objective-function selector for the Analysis tab.

Lets the user steer what the simulation optimizes, in plain language:

* **Maximize growth** — the biomass reaction (auto-detected).
* **Maximize a product** — pick any reaction (e.g. a product exchange).
* **Balance growth & product** — maximize ``growth + weight × product`` so the
  user can trade biomass against a target chemical.

Emits ``objective_changed(terms, direction)`` where ``terms`` is
``{reaction_id: weight}`` for the core to apply as the model objective.
"""

from __future__ import annotations

from typing import List, Optional

import cobra
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ...core import editing

MODE_GROWTH = "growth"
MODE_PRODUCT = "product"
MODE_BALANCE = "balance"


class ObjectiveBar(QGroupBox):
    objective_changed = Signal(dict, str)  # {reaction_id: weight}, direction
    consortia_requested = Signal()         # open the consortia objective settings

    def __init__(self):
        super().__init__("Optimization objective")
        self._biomass: Optional[str] = None
        self._model: Optional[cobra.Model] = None

        self.current_label = QLabel("Current objective: —")
        self.current_label.setWordWrap(True)
        self.current_label.setStyleSheet("color: #6B7280;")

        self.mode_combo = QComboBox()
        self.mode_combo.addItem("Maximize growth (biomass)", MODE_GROWTH)
        self.mode_combo.addItem("Maximize a product", MODE_PRODUCT)
        self.mode_combo.addItem("Balance growth & product", MODE_BALANCE)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)

        self.product_combo = QComboBox()
        self.product_combo.setEditable(True)
        self.product_combo.setMinimumWidth(180)

        self.weight = QDoubleSpinBox()
        self.weight.setRange(0.0, 1000.0)
        self.weight.setDecimals(2)
        self.weight.setSingleStep(0.1)
        self.weight.setValue(0.5)
        self.weight.setToolTip("Relative reward for the product vs. growth (growth weight = 1).")

        self.apply_btn = QPushButton("Apply objective")
        self.apply_btn.setObjectName("primary")
        self.apply_btn.clicked.connect(self._emit)

        self.custom_btn = QPushButton("Custom / weighted objective…")
        self.custom_btn.setToolTip(
            "Build an objective from several reactions with weights "
            "(e.g. 60% biomass + 40% product).")
        self.custom_btn.clicked.connect(self._open_custom)

        # Shown only for community/consortium models: dominance + no-starvation floor.
        self.consortia_btn = QPushButton("Consortia objective (dominance & min growth)…")
        self.consortia_btn.setToolTip(
            "This is a community model. Set each organism's dominance and a minimum "
            "growth so the faster grower can't starve the others.")
        self.consortia_btn.clicked.connect(self.consortia_requested)
        self.consortia_btn.setVisible(False)

        self.product_label = QLabel("Product:")
        self.weight_label = QLabel("Product weight:")

        row = QHBoxLayout()
        row.addWidget(self.product_label)
        row.addWidget(self.product_combo, 1)
        row.addWidget(self.weight_label)
        row.addWidget(self.weight)

        layout = QVBoxLayout(self)
        layout.addWidget(self.current_label)
        layout.addWidget(self.mode_combo)
        layout.addLayout(row)
        layout.addWidget(self.apply_btn)
        layout.addWidget(self.custom_btn)
        layout.addWidget(self.consortia_btn)

        self._on_mode_changed()

    def set_community(self, active: bool) -> None:
        """Toggle the consortia-objective button (shown only for community models)."""
        self.consortia_btn.setVisible(bool(active))

    def set_model(self, model: Optional[cobra.Model]) -> None:
        self._model = model
        self.product_combo.clear()
        if model is None:
            self.current_label.setText("Current objective: —")
            return
        self._biomass = editing.guess_biomass_reaction(model)
        ids = [r.id for r in model.reactions]
        self.product_combo.addItems(ids)
        # Default the product to an exchange reaction if present.
        for rid in ids:
            if rid.startswith("EX_"):
                self.product_combo.setCurrentText(rid)
                break
        self.update_current(model)

    def update_current(self, model: cobra.Model) -> None:
        terms = editing.current_objective_terms(model)
        if not terms:
            self.current_label.setText("Current objective: (none set)")
            return
        pretty = "  +  ".join(f"{w:g}·{rid}" for rid, w in terms.items())
        direction = getattr(model, "objective_direction", "max")
        self.current_label.setText(f"Current objective: {direction}  {pretty}")

    def _on_mode_changed(self) -> None:
        mode = self.mode_combo.currentData()
        show_product = mode in (MODE_PRODUCT, MODE_BALANCE)
        for w in (self.product_label, self.product_combo):
            w.setVisible(show_product)
        show_weight = mode == MODE_BALANCE
        for w in (self.weight_label, self.weight):
            w.setVisible(show_weight)

    def _open_custom(self) -> None:
        if self._model is None:
            return
        from ..dialogs.objective_dialog import ObjectiveDialog

        dlg = ObjectiveDialog(self._model, self)
        if dlg.exec() == ObjectiveDialog.Accepted:
            coeffs, direction = dlg.result()
            if coeffs:
                self.objective_changed.emit(coeffs, direction)

    def _emit(self) -> None:
        mode = self.mode_combo.currentData()
        product = self.product_combo.currentText().strip()
        if mode == MODE_GROWTH:
            if not self._biomass:
                return
            terms = {self._biomass: 1.0}
        elif mode == MODE_PRODUCT:
            if not product:
                return
            terms = {product: 1.0}
        else:  # balance
            terms = {}
            if self._biomass:
                terms[self._biomass] = 1.0
            if product:
                terms[product] = self.weight.value()
        if terms:
            self.objective_changed.emit(terms, "max")
