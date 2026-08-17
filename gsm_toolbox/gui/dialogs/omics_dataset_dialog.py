"""Dataset Preparation wizard: map a raw omics table onto the model's namespace.

Handles transcriptomics / proteomics (mapped to genes) and metabolomics (mapped to
metabolites). The user picks the file, the data type, the identifier column and the
value column(s) (replicate samples can be averaged), then maps onto the model and
sees a coverage summary before loading the result.
"""

from __future__ import annotations

from typing import Optional

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)

from ...core import omics_prep as op

_KINDS = [
    ("transcriptomics", "Transcriptomics (RNA-seq / microarray → genes)"),
    ("proteomics", "Proteomics (protein abundance → genes)"),
    ("metabolomics", "Metabolomics (metabolite levels → metabolites)"),
]
_AGGS = ["mean", "median", "max", "min", "sum", "first"]


class OmicsDatasetDialog(QDialog):
    def __init__(self, parent, model: cobra.Model):
        super().__init__(parent)
        self.setWindowTitle("Prepare omics dataset")
        self.resize(680, 640)
        self._model = model
        self._df = None
        self.dataset: Optional[op.PreparedDataset] = None

        v = QVBoxLayout(self)
        v.addWidget(QLabel(
            "Load a raw omics table in any common format (CSV / TSV / Excel). The tool maps "
            "its identifiers onto this model and produces the (id, value) data used for "
            "context-specific analysis."))

        # --- file + kind ---
        fbox = QGroupBox("1 · Source file and data type")
        fv = QVBoxLayout(fbox)
        row = QHBoxLayout()
        self.browse_btn = QPushButton("Choose file…")
        self.browse_btn.clicked.connect(self._choose_file)
        self.path_label = QLabel("<i>no file chosen</i>")
        self.path_label.setWordWrap(True)
        row.addWidget(self.browse_btn)
        row.addWidget(self.path_label, 1)
        fv.addLayout(row)
        krow = QHBoxLayout()
        krow.addWidget(QLabel("Data type:"))
        self.kind_combo = QComboBox()
        for k, lbl in _KINDS:
            self.kind_combo.addItem(lbl, k)
        krow.addWidget(self.kind_combo, 1)
        fv.addLayout(krow)
        v.addWidget(fbox)

        # --- columns ---
        self.cbox = QGroupBox("2 · Columns")
        cv = QVBoxLayout(self.cbox)
        idrow = QHBoxLayout()
        idrow.addWidget(QLabel("Identifier column:"))
        self.id_combo = QComboBox()
        idrow.addWidget(self.id_combo, 1)
        idrow.addWidget(QLabel("Combine samples by:"))
        self.agg_combo = QComboBox()
        self.agg_combo.addItems(_AGGS)
        idrow.addWidget(self.agg_combo)
        cv.addLayout(idrow)
        cv.addWidget(QLabel("Value column(s) — tick the samples/replicates to use:"))
        self.value_list = QListWidget()
        self.value_list.setMaximumHeight(140)
        cv.addWidget(self.value_list)
        self.map_btn = QPushButton("Map to model")
        self.map_btn.setObjectName("primary")
        self.map_btn.clicked.connect(self._map)
        cv.addWidget(self.map_btn)
        self.cbox.setEnabled(False)
        v.addWidget(self.cbox)

        # --- summary + preview ---
        sbox = QGroupBox("3 · Mapping summary")
        sv = QVBoxLayout(sbox)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setPlaceholderText("Map to the model to see how well the identifiers matched.")
        self.summary.setMaximumHeight(200)
        sv.addWidget(self.summary)
        v.addWidget(sbox)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.button(QDialogButtonBox.Ok).setText("Use this dataset")
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(False)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        # extra: save the prepared table to disk
        self.save_btn = self.buttons.addButton("Save prepared table…", QDialogButtonBox.ActionRole)
        self.save_btn.setEnabled(False)
        self.save_btn.clicked.connect(self._save)
        v.addWidget(self.buttons)

    # ---- steps ----
    def _choose_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose omics table", "",
            "Omics tables (*.csv *.tsv *.txt *.xlsx *.xls);;All files (*)")
        if not path:
            return
        try:
            self._df = op.read_table(path)
        except op.OmicsPrepError as exc:
            QMessageBox.critical(self, "Could not read file", str(exc))
            return
        self.path_label.setText(path)
        self._populate_columns()
        self.cbox.setEnabled(True)
        self.summary.clear()
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(False)
        self.save_btn.setEnabled(False)

    def _populate_columns(self) -> None:
        cols = [str(c) for c in self._df.columns]
        numeric = set(op.numeric_columns(self._df))
        self.id_combo.clear()
        self.id_combo.addItems(cols)
        self.id_combo.setCurrentText(op.guess_id_column(self._df))
        self.value_list.clear()
        for c in cols:
            it = QListWidgetItem(c + ("  (numeric)" if c in numeric else ""))
            it.setData(Qt.UserRole, c)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if c in numeric and c != self.id_combo.currentText()
                             else Qt.Unchecked)
            self.value_list.addItem(it)

    def _checked_values(self) -> list:
        return [self.value_list.item(i).data(Qt.UserRole)
                for i in range(self.value_list.count())
                if self.value_list.item(i).checkState() == Qt.Checked]

    def _map(self) -> None:
        if self._df is None:
            return
        value_cols = self._checked_values()
        if not value_cols:
            QMessageBox.information(self, "Pick value columns",
                                    "Tick at least one numeric value/sample column.")
            return
        try:
            ds = op.prepare_dataset(
                self._df, self._model, kind=self.kind_combo.currentData(),
                id_column=self.id_combo.currentText(), value_columns=value_cols,
                aggregate=self.agg_combo.currentText())
        except op.OmicsPrepError as exc:
            QMessageBox.critical(self, "Could not map dataset", str(exc))
            return
        self.dataset = ds
        text = ds.summary.text()
        if ds.summary.n_model_targets == 0:
            text += ("\n\n⚠ Nothing matched. Check that the identifier column and data type "
                     "are right, and that the ids use the same scheme as the model's genes.")
        elif ds.summary.coverage < 0.15:
            text += ("\n\n⚠ Low coverage — the identifier scheme may differ from the model's. "
                     "Try a different id column, or a source with matching gene ids.")
        self.summary.setPlainText(text)
        ok = ds.summary.n_model_targets > 0
        self.buttons.button(QDialogButtonBox.Ok).setEnabled(ok)
        self.save_btn.setEnabled(ok)

    def _save(self) -> None:
        if self.dataset is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save prepared table", "prepared_dataset.csv",
            "CSV (*.csv);;TSV (*.tsv)")
        if not path:
            return
        try:
            op.write_prepared_csv(self.dataset, path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Could not save", str(exc))
            return
        QMessageBox.information(self, "Saved", f"Prepared table written to:\n{path}")
