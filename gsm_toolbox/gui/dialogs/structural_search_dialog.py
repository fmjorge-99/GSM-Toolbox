"""Structural Search — find database compounds related to a target you describe (VI.2).

The target box matches names literally, which is fast and predictable but hides compounds
stored under a different name: Δ⁹-THC is absent while its family sits under *Cannabidiolic
acid*; violacein is stored as *violaceinate*. Rather than making the main search fuzzy —
which would trade one surprise for another — this dialog is an explicit second route in:
give a **name**, a **SMILES** or an **InChIKey**, and it lists the compounds in the loaded
databases that are structurally or chemically related, with the number of reactions that
can actually produce each one.
"""
from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox, QDialog, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout)

from .. import style
from ..widgets.busy import run_busy, was_cancelled


class StructuralSearchDialog(QDialog):
    """Emits :attr:`target_chosen` with a metabolite id when the user picks a row."""

    target_chosen = Signal(str)

    def __init__(self, parent, database, *, host_name: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Structural search")
        self.resize(820, 560)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
        self._db = database
        self._hits = []

        v = QVBoxLayout(self)
        v.addWidget(QLabel(
            "<b>Find a target by structure or chemistry.</b><br>"
            "Use this when the exact name is not in the database — many compounds are "
            "stored as an acid, an anion or a homologue."))

        row = QHBoxLayout()
        self.kind = QComboBox()
        self.kind.addItems(["Compound name", "SMILES", "InChIKey"])
        row.addWidget(self.kind)
        self.query = QLineEdit()
        self.query.setPlaceholderText("Compound name, SMILES or InChIKey")
        self.query.returnPressed.connect(self._search)
        row.addWidget(self.query, 1)
        self.go = QPushButton("Search")
        self.go.setObjectName("primary")
        self.go.clicked.connect(self._search)
        row.addWidget(self.go)
        v.addLayout(row)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setStyleSheet(f"color:{style.TEXT_MUTED};")
        v.addWidget(self.status)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["Compound", "Metabolite id", "Producing reactions", "Why it matched"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.itemDoubleClicked.connect(lambda *_: self._use())
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.Stretch)
        v.addWidget(self.table, 1)

        v.addWidget(QLabel(
            "<span style='color:#5f6368'>A compound with <b>0 producing reactions</b> is "
            "present but nothing can make it — prefer a candidate with several.</span>"))

        btns = QHBoxLayout()
        self.use_btn = QPushButton("Use as target")
        self.use_btn.setObjectName("primary")
        self.use_btn.setEnabled(False)
        self.use_btn.clicked.connect(self._use)
        btns.addWidget(self.use_btn)
        btns.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        btns.addWidget(close)
        v.addLayout(btns)

    # -- search ---------------------------------------------------------------------
    def _search(self) -> None:
        q = self.query.text().strip()
        if not q or self._db is None:
            return
        kind = self.kind.currentText()

        def work():
            from ...core import pathway_search as ps
            name, smiles = q, ""
            if kind == "SMILES":
                smiles, name = q, _name_for_smiles(q)
            elif kind == "InChIKey":
                smiles, name = _smiles_for_inchikey(q), _name_for_inchikey(q) or q
            else:
                smiles = _smiles_for_name(q)
            return ps.nearest_reachable_analogues(self._db, name, target_smiles=smiles,
                                                  limit=25), name, smiles

        ok, res = run_busy(self, f"Searching the loaded databases for “{q}”…", work,
                           title="Structural search", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                self.status.setText(f"Search failed: {res}")
            return
        hits, name, smiles = res
        self._hits = hits
        self.table.setRowCount(0)
        for h in hits:
            r = self.table.rowCount()
            self.table.insertRow(r)
            self.table.setItem(r, 0, QTableWidgetItem(h["name"]))
            self.table.setItem(r, 1, QTableWidgetItem(h["id"]))
            item = QTableWidgetItem(str(h["n_producers"]))
            item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(r, 2, item)
            self.table.setItem(r, 3, QTableWidgetItem(h["reason"]))
        self.use_btn.setEnabled(bool(hits))
        if hits:
            self.table.selectRow(0)
            extra = f" (structure resolved: {smiles[:40]})" if smiles else ""
            self.status.setText(f"{len(hits)} related compound(s) found for "
                                f"“{name}”{extra}. Double-click one to use it.")
        else:
            self.status.setText(
                f"No related compound found for “{q}”. Try a shorter fragment of the "
                "name, or a different identifier — and check the right databases are "
                "loaded.")

    def _use(self) -> None:
        row = self.table.currentRow()
        if 0 <= row < len(self._hits):
            self.target_chosen.emit(self._hits[row]["id"])
            self.accept()


# -- identifier → structure helpers (all tolerant of being offline) --------------------
def _smiles_for_name(name: str) -> str:
    try:
        from ..widgets.structure_fetcher import _pubchem_smiles, _ambiguous_name
        if name and not _ambiguous_name(name):
            return _pubchem_smiles(name=name) or ""
    except Exception:  # noqa: BLE001
        pass
    return ""


def _smiles_for_inchikey(ik: str) -> str:
    try:
        from ..widgets.structure_fetcher import _pubchem_smiles
        return _pubchem_smiles(inchikey=ik.strip()) or ""
    except Exception:  # noqa: BLE001
        return ""


def _name_for_inchikey(ik: str) -> str:
    try:
        from ..widgets.structure_fetcher import name_from_inchikey
        return name_from_inchikey(ik.strip()) or ""
    except Exception:  # noqa: BLE001
        return ""


def _name_for_smiles(smiles: str) -> str:
    try:
        from ..widgets.structure_fetcher import name_from_smiles
        return name_from_smiles(smiles) or smiles
    except Exception:  # noqa: BLE001
        return smiles
