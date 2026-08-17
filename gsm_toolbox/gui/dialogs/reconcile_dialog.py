"""Review the identifier correspondences the toolbox proposes between model and database.

A wrong automatic merge is worse than no merge. Two metabolites fused by mistake change
what the model can make, and nothing downstream looks wrong: the route balances, carries
flux and reports a yield. So the software proposes and the user decides.

Every row states the score and the evidence behind it in the same words the core module
produced, and a pair the checks rejected outright is shown greyed with its reason rather
than hidden. A rejected pair the user recognises as correct is exactly what an override
is for, and the override is recorded as one.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QDialog, QHBoxLayout, QHeaderView, QLabel,
    QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from ... core import id_reconcile as rec
from .. import style

_COLUMNS = ["Use", "Database metabolite", "Treated as", "Confidence", "Score", "Why"]


class ReconcileDialog(QDialog):
    """Accept or reject proposed metabolite correspondences."""

    def __init__(self, parent, proposals: Sequence[rec.Proposal], *,
                 database_name: str = "", preselected: Optional[Dict[str, str]] = None):
        super().__init__(parent)
        self.setWindowTitle("Reconcile metabolite identifiers")
        self.resize(1080, 620)
        self._proposals = list(proposals)
        self._overrides: Dict[str, str] = {}

        outer = QVBoxLayout(self)

        heading = QLabel(
            "The model and the database can use different identifiers for the same "
            "compound. Where that happens the database metabolite is added as an "
            "orphan and the route carries no flux, even though it looks complete.\n\n"
            "Tick the correspondences you want applied. High confidence rows are "
            "ticked already. Rows the chemistry checks rejected are shown in grey and "
            "need a deliberate override.")
        heading.setWordWrap(True)
        outer.addWidget(heading)

        self.table = QTableWidget(len(self._proposals), len(_COLUMNS))
        self.table.setHorizontalHeaderLabels(_COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(len(_COLUMNS) - 1, QHeaderView.Stretch)

        self._boxes: List[QCheckBox] = []
        for row, proposal in enumerate(self._proposals):
            box = QCheckBox()
            chosen = (preselected or {}).get(proposal.db_id) == proposal.host_id
            box.setChecked(bool(chosen) or proposal.auto_acceptable)
            box.setEnabled(not proposal.blocked)
            if proposal.blocked:
                box.setToolTip("Blocked: " + proposal.blocked +
                               "\nUse Override to apply it anyway.")
            self._boxes.append(box)
            holder = QWidget()
            layout = QHBoxLayout(holder)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.addWidget(box)
            layout.setAlignment(Qt.AlignCenter)
            self.table.setCellWidget(row, 0, holder)

            values = [
                f"{proposal.db_id}"
                + (f"  ({proposal.db_name})" if proposal.db_name else ""),
                f"{proposal.host_id}"
                + (f"  ({proposal.host_name})" if proposal.host_name else ""),
                proposal.confidence,
                f"{proposal.score:.0f}",
                proposal.why(),
            ]
            for column, text in enumerate(values, start=1):
                item = QTableWidgetItem(text)
                item.setToolTip(text)
                if proposal.blocked:
                    item.setForeground(Qt.gray)
                self.table.setItem(row, column, item)
        self.table.resizeColumnsToContents()
        outer.addWidget(self.table, 1)

        self.status = QLabel("")
        self.status.setStyleSheet(f"color:{style.TEXT_MUTED};")
        outer.addWidget(self.status)

        buttons = QHBoxLayout()
        select_high = QPushButton("Select high confidence")
        select_high.clicked.connect(self._select_high)
        buttons.addWidget(select_high)
        clear = QPushButton("Clear all")
        clear.clicked.connect(lambda: self._set_all(False))
        buttons.addWidget(clear)
        override = QPushButton("Override selected row")
        override.setToolTip("Apply a correspondence the chemistry checks rejected. "
                            "The reason is recorded with the mapping.")
        override.clicked.connect(self._override_selected)
        buttons.addWidget(override)
        buttons.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(cancel)
        apply_button = QPushButton("Apply mapping")
        apply_button.setObjectName("primary")
        apply_button.clicked.connect(self.accept)
        buttons.addWidget(apply_button)
        outer.addLayout(buttons)

        self._refresh_status()
        for box in self._boxes:
            box.toggled.connect(self._refresh_status)

        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

    # -- helpers -----------------------------------------------------------------------
    def _set_all(self, state: bool) -> None:
        for box, proposal in zip(self._boxes, self._proposals):
            if not proposal.blocked:
                box.setChecked(state)

    def _select_high(self) -> None:
        for box, proposal in zip(self._boxes, self._proposals):
            if not proposal.blocked:
                box.setChecked(proposal.auto_acceptable)

    def _override_selected(self) -> None:
        row = self.table.currentRow()
        if row < 0:
            return
        proposal = self._proposals[row]
        if not proposal.blocked:
            QMessageBox.information(self, "Override",
                                    "This row is not blocked, so it needs no override.")
            return
        answer = QMessageBox.question(
            self, "Override a rejected correspondence",
            f"{proposal.db_id} was rejected because {proposal.blocked}.\n\n"
            f"{proposal.why()}\n\n"
            f"Applying it anyway will treat these as the same compound everywhere. "
            f"Record the override and continue?",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel)
        if answer != QMessageBox.Yes:
            return
        self._overrides[proposal.db_id] = proposal.blocked
        self._boxes[row].setEnabled(True)
        self._boxes[row].setChecked(True)
        for column in range(1, len(_COLUMNS)):
            item = self.table.item(row, column)
            if item is not None:
                item.setForeground(Qt.black)
        self._refresh_status()

    def _refresh_status(self) -> None:
        chosen = sum(1 for b in self._boxes if b.isChecked())
        blocked = sum(1 for p in self._proposals if p.blocked)
        text = f"{chosen} of {len(self._proposals)} correspondences selected"
        if blocked:
            text += f", {blocked} rejected by the chemistry checks"
        if self._overrides:
            text += f", {len(self._overrides)} overridden"
        self.status.setText(text)

    # -- results -----------------------------------------------------------------------
    def accepted(self) -> List[rec.Proposal]:
        return [p for p, b in zip(self._proposals, self._boxes) if b.isChecked()]

    def mapping(self) -> Dict[str, str]:
        return {p.db_id: p.host_id for p in self.accepted()}

    def record(self, database_name: str = "") -> dict:
        """The versioned artefact to store in the project."""
        return rec.to_record(self.accepted(), database=database_name,
                             overrides=self._overrides)
