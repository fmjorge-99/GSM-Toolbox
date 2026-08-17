"""The Results tab: a header summary plus a sortable, exportable table.

Generic enough to show any pandas DataFrame (FBA fluxes, FVA ranges, deletion
screens). Phase 1 uses it for FBA/pFBA flux tables.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableView,
    QVBoxLayout,
    QWidget,
)


class _DataFrameModel(QAbstractTableModel):
    def __init__(self, df: Optional[pd.DataFrame] = None):
        super().__init__()
        self._df = df if df is not None else pd.DataFrame()
        self._editable_cols: set = set()   # column NAMES the user may edit (#B5)

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        # Keep the SAME object (edits must persist back to the caller's DataFrame).
        self._df = df
        self.endResetModel()

    def set_editable_columns(self, columns) -> None:
        self._editable_cols = set(columns or ())

    def dataframe(self) -> pd.DataFrame:
        return self._df

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._df)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return len(self._df.columns)

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return str(self._df.columns[section])
        return str(self._df.index[section])

    def _col_editable(self, col: int) -> bool:
        return (0 <= col < len(self._df.columns)
                and str(self._df.columns[col]) in self._editable_cols)

    def flags(self, index):  # noqa: N802
        base = super().flags(index)
        if index.isValid() and self._col_editable(index.column()):
            return base | Qt.ItemIsEditable
        return base

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or role not in (Qt.DisplayRole, Qt.EditRole):
            return None
        value = self._df.iat[index.row(), index.column()]
        if role == Qt.EditRole:
            return str(value)
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)

    def setData(self, index, value, role=Qt.EditRole):  # noqa: N802
        if not index.isValid() or role != Qt.EditRole or not self._col_editable(index.column()):
            return False
        text = str(value).strip()
        if not text:
            return False
        self._df.iat[index.row(), index.column()] = text
        self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.EditRole])
        return True


class ResultsView(QWidget):
    # Emitted on right-click of a row: the row's cell strings + global position.
    # The main window scans the cells for a reaction id and builds an actions menu.
    row_context_requested = Signal(list, object)

    def __init__(self):
        super().__init__()
        self.header = QLabel("No analysis run yet.")
        self.header.setWordWrap(True)
        self.header.setStyleSheet("font-weight: 600; padding: 4px;")

        self._model = _DataFrameModel()
        self.table = QTableView()
        self.table.setModel(self._model)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._emit_row_context)

        self.export_btn = QPushButton("Export to CSV…")
        self.export_btn.clicked.connect(self._export)
        self.export_btn.setEnabled(False)

        top = QHBoxLayout()
        top.addWidget(self.header, 1)
        top.addWidget(self.export_btn)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.table, 1)

        self._df: Optional[pd.DataFrame] = None

    def show_dataframe(self, df: pd.DataFrame, header: str) -> None:
        self._df = df
        self._header = header
        self._model.set_dataframe(df)
        self.header.setText(header)
        self.table.resizeColumnsToContents()
        self.export_btn.setEnabled(not df.empty)

    def set_editable_columns(self, columns) -> None:
        """Make the named columns user-editable (e.g. 'suggested_id') (#B5)."""
        from PySide6.QtWidgets import QAbstractItemView
        self._model.set_editable_columns(columns)
        self.table.setEditTriggers(QAbstractItemView.DoubleClicked
                                   | QAbstractItemView.EditKeyPressed)

    def current_dataframe(self) -> Optional[pd.DataFrame]:
        return self._df

    def current_header(self) -> str:
        return getattr(self, "_header", "")

    def _emit_row_context(self, pos) -> None:
        index = self.table.indexAt(pos)
        if not index.isValid() or self._df is None or self._df.empty:
            return
        row = index.row()
        try:
            cells = [str(v) for v in self._df.iloc[row].tolist()]
        except Exception:  # noqa: BLE001
            return
        self.row_context_requested.emit(cells, self.table.viewport().mapToGlobal(pos))

    def _export(self) -> None:
        if self._df is None:
            return
        from ..widgets.dialog_util import choose_save_path
        path = choose_save_path(self, "Export results", "results.csv", "CSV (*.csv)")
        if path:
            if not path.lower().endswith(".csv"):
                path += ".csv"
            self._df.to_csv(path, index=False)
