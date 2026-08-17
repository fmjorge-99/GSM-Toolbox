"""Growth medium / constraints editor.

Lists every exchange reaction with editable uptake (max import) bounds, plus an
aerobic/anaerobic toggle. Returns the new medium as ``{exchange_id: max_uptake}``.
"""

from __future__ import annotations

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from ...core import media


class MediaDialog(QDialog):
    def __init__(self, model: cobra.Model, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Growth Medium / Constraints")
        self.setMinimumSize(560, 480)
        self._model = model

        info = QLabel(
            "Set the maximum uptake rate (mmol/gDW/h) for each nutrient the cell can "
            "import. A value of 0 removes the nutrient from the medium."
        )
        info.setWordWrap(True)

        self.aerobic = QCheckBox("Aerobic (oxygen available)")
        o2 = media.find_oxygen_exchange(model)
        if o2 is not None:
            self.aerobic.setChecked(model.reactions.get_by_id(o2).lower_bound < 0)
        else:
            self.aerobic.setEnabled(False)
            self.aerobic.setToolTip("No oxygen exchange reaction found in this model.")

        rows = media.list_exchanges(model)
        self.table = QTableWidget(len(rows), 3)
        self.table.setHorizontalHeaderLabels(["Exchange", "Name", "Max uptake"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
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

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(info)
        layout.addWidget(self.aerobic)
        layout.addWidget(self.table, 1)
        layout.addWidget(buttons)

    def medium(self) -> dict:
        """Return only exchanges with a positive uptake (the active medium)."""
        return {rid: spin.value() for rid, spin in self._spins.items() if spin.value() > 0}

    def is_aerobic(self) -> bool:
        return self.aerobic.isChecked()


def _ro_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item
