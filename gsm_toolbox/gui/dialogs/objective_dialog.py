"""Custom multi-term optimization objective editor.

Lets the user build an objective from any number of reactions, each with a weight,
e.g. 60% biomass + 40% product. Weights can be treated as *relative priorities*
(normalized by each reaction's maximum attainable flux) so the percentages behave
intuitively regardless of the reactions' raw flux magnitudes.

Returns objective coefficients ``{reaction_id: coefficient}`` and a direction.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from ...core import editing


class ObjectiveDialog(QDialog):
    def __init__(self, model: cobra.Model, parent=None):
        super().__init__(parent)
        self._model = model
        self._reaction_ids = [r.id for r in model.reactions]
        self.setWindowTitle("Custom optimization objective")
        self.setMinimumSize(560, 420)

        intro = QLabel(
            "Build an objective from one or more reactions, each with a weight. "
            "Example: 0.6 for the biomass reaction and 0.4 for a product exchange to "
            "balance growth with production.")
        intro.setWordWrap(True)

        self.direction = QComboBox()
        self.direction.addItem("Maximize", "max")
        self.direction.addItem("Minimize", "min")
        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Direction:"))
        dir_row.addWidget(self.direction)
        dir_row.addStretch(1)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Reaction", "Weight"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.verticalHeader().setVisible(False)

        add_btn = QPushButton("Add reaction")
        add_btn.clicked.connect(lambda: self._add_row())
        rm_btn = QPushButton("Remove selected")
        rm_btn.clicked.connect(self._remove_row)
        btn_row = QHBoxLayout()
        btn_row.addWidget(add_btn)
        btn_row.addWidget(rm_btn)
        btn_row.addStretch(1)

        self.normalize = QCheckBox(
            "Treat weights as relative priorities (normalize by each reaction's max flux)")
        self.normalize.setChecked(True)
        self.normalize.setToolTip(
            "Recommended when mixing reactions with very different flux scales "
            "(e.g. growth ≈ 1 vs a product exchange ≈ 20).")

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(dir_row)
        layout.addWidget(self.table, 1)
        layout.addLayout(btn_row)
        layout.addWidget(self.normalize)
        layout.addWidget(buttons)

        self._prefill()

    def _prefill(self) -> None:
        current = editing.current_objective_terms(self._model)
        if not current:
            biomass = editing.guess_biomass_reaction(self._model)
            current = {biomass: 1.0} if biomass else {}
        for rid, coeff in current.items():
            self._add_row(rid, abs(coeff) if coeff else 1.0)
        if not current:
            self._add_row()

    def _add_row(self, reaction_id: Optional[str] = None, weight: float = 1.0) -> None:
        row = self.table.rowCount()
        self.table.insertRow(row)
        combo = QComboBox()
        combo.setEditable(True)
        combo.addItems(self._reaction_ids)
        if reaction_id and reaction_id in self._reaction_ids:
            combo.setCurrentText(reaction_id)
        self.table.setCellWidget(row, 0, combo)
        spin = QDoubleSpinBox()
        spin.setRange(0.0, 1e6)
        spin.setDecimals(3)
        spin.setSingleStep(0.1)
        spin.setValue(weight)
        self.table.setCellWidget(row, 1, spin)

    def _remove_row(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.table.removeRow(r)

    def _raw_weights(self) -> Dict[str, float]:
        weights: Dict[str, float] = {}
        for row in range(self.table.rowCount()):
            combo = self.table.cellWidget(row, 0)
            spin = self.table.cellWidget(row, 1)
            rid = combo.currentText().strip()
            w = spin.value()
            if rid and w:
                weights[rid] = weights.get(rid, 0.0) + w
        return weights

    def result(self) -> Tuple[Dict[str, float], str]:
        """Return (coefficients, direction). Coefficients are normalized if requested."""
        weights = self._raw_weights()
        if self.normalize.isChecked() and len(weights) > 1:
            coeffs = editing.normalized_objective_coefficients(self._model, weights)
        else:
            coeffs = weights
        return coeffs, self.direction.currentData()
