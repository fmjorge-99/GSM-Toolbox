"""The Omics tab: shows the omics datasets mapped onto the model.

This tab is added to the main window only once an omics dataset has been prepared
or loaded (a mapping to the model). It lists each dataset, its coverage summary,
and the (id, value) values, so the user can see exactly what feeds eFlux / GIMME.
"""

from __future__ import annotations

from typing import Dict

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)


class OmicsPanel(QWidget):
    def __init__(self):
        super().__init__()
        self._datasets: Dict[str, dict] = {}   # name -> {values, summary, kind}

        top = QHBoxLayout()
        top.addWidget(QLabel("<b>Omics datasets mapped to the model</b>"))
        top.addStretch(1)
        top.addWidget(QLabel("Dataset:"))
        self.selector = QComboBox()
        self.selector.setMinimumWidth(240)
        self.selector.currentIndexChanged.connect(self._show_current)
        top.addWidget(self.selector)

        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setMaximumHeight(150)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["id", "value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSortingEnabled(True)

        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.summary)
        lay.addWidget(self.table, 1)

    def add_dataset(self, name: str, values: Dict[str, float], *, summary: str = "",
                    kind: str = "") -> None:
        """Add (or replace) a dataset and show it."""
        self._datasets[name] = {"values": dict(values), "summary": summary, "kind": kind}
        self.selector.blockSignals(True)
        if self.selector.findText(name) < 0:
            self.selector.addItem(name)
        self.selector.setCurrentText(name)
        self.selector.blockSignals(False)
        self._show_current()

    def _show_current(self) -> None:
        name = self.selector.currentText()
        data = self._datasets.get(name)
        if not data:
            return
        kind = data.get("kind", "")
        id_header = "metabolite" if kind == "metabolomics" else ("gene" if kind else "id")
        self.summary.setPlainText(data.get("summary") or
                                  f"{len(data['values'])} values loaded.")
        values = data["values"]
        self.table.setSortingEnabled(False)
        self.table.setHorizontalHeaderLabels([id_header, "value"])
        self.table.setRowCount(len(values))
        for row, (k, v) in enumerate(sorted(values.items())):
            self.table.setItem(row, 0, QTableWidgetItem(str(k)))
            vi = QTableWidgetItem()
            try:
                vi.setData(Qt.EditRole, float(v))
            except (TypeError, ValueError):
                vi.setText(str(v))
            self.table.setItem(row, 1, vi)
        self.table.setSortingEnabled(True)

    def clear_all(self) -> None:
        self._datasets.clear()
        self.selector.blockSignals(True)
        self.selector.clear()
        self.selector.blockSignals(False)
        self.summary.clear()
        self.table.setRowCount(0)
