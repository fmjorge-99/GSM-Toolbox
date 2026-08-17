"""Dual-panel substrate picker (Issue R4).

A quick shortcut for choosing which nutrients/substrates the cell may take up:
the left panel lists every exchangeable metabolite in the model (searchable), the
right panel holds the currently-selected uptakes. "Add >" / "< Remove" move items
between the two. Finer control (per-substrate uptake rate) still lives in the
Growth Settings "Medium" table.

Returns the chosen set of exchange ids via :meth:`selected_ids`.
"""

from __future__ import annotations

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
)

from ...core import physiology


class SubstratePickerDialog(QDialog):
    def __init__(self, model: cobra.Model, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select growth substrates")
        self.setMinimumSize(720, 480)
        self._rows = physiology.substrate_exchanges(model)
        self._by_id = {r["id"]: r for r in self._rows}

        info = QLabel("Choose which nutrients the cell may take up. Move carbon sources "
                      "and other substrates into “Selected”. Uptake rates can be tuned "
                      "afterwards in the Medium table.")
        info.setWordWrap(True)

        # Left: available (searchable)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search metabolites…")
        self.search.textChanged.connect(self._filter)
        self.available = QListWidget()
        self.available.setSelectionMode(QAbstractItemView.ExtendedSelection)
        left = QVBoxLayout()
        left.addWidget(QLabel("Available substrates"))
        left.addWidget(self.search)
        left.addWidget(self.available, 1)

        # Middle: add/remove
        self.add_btn = QPushButton("Add >")
        self.remove_btn = QPushButton("< Remove")
        self.add_btn.clicked.connect(self._add)
        self.remove_btn.clicked.connect(self._remove)
        mid = QVBoxLayout()
        mid.addStretch(1)
        mid.addWidget(self.add_btn)
        mid.addWidget(self.remove_btn)
        mid.addStretch(1)

        # Right: selected
        self.selected = QListWidget()
        self.selected.setSelectionMode(QAbstractItemView.ExtendedSelection)
        right = QVBoxLayout()
        right.addWidget(QLabel("Selected for uptake"))
        right.addWidget(self.selected, 1)

        cols = QHBoxLayout()
        cols.addLayout(left, 5)
        cols.addLayout(mid, 1)
        cols.addLayout(right, 5)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(info)
        layout.addLayout(cols, 1)
        layout.addWidget(buttons)

        # Seed the two panels from current uptake state.
        for r in self._rows:
            (self.selected if r["active"] else self.available).addItem(self._make_item(r))
        self._filter(self.search.text())

    def _make_item(self, row: dict) -> QListWidgetItem:
        label = f"{row['name']}  ({row['id']})"
        if row["carbon"]:
            label += f"  · C{row['carbon']}"
        item = QListWidgetItem(label)
        item.setData(Qt.UserRole, row["id"])
        return item

    def _filter(self, text: str) -> None:
        text = (text or "").lower()
        for i in range(self.available.count()):
            it = self.available.item(i)
            it.setHidden(bool(text) and text not in it.text().lower())

    def _move(self, src: QListWidget, dst: QListWidget) -> None:
        for it in src.selectedItems():
            rid = it.data(Qt.UserRole)
            src.takeItem(src.row(it))
            dst.addItem(self._make_item(self._by_id[rid]))
        self._filter(self.search.text())

    def _add(self) -> None:
        self._move(self.available, self.selected)

    def _remove(self) -> None:
        self._move(self.selected, self.available)

    def selected_ids(self) -> list:
        return [self.selected.item(i).data(Qt.UserRole) for i in range(self.selected.count())]
