"""Suggest EC numbers for a reaction, and the enzyme sequences that could perform it.

Most reactions in a universal database carry no EC annotation, which leaves the user with
no way in to "which enzyme do I express?". This dialog answers that in two stages:

1. **EC numbers** — annotated, or looked up from the reaction's own KEGG/Rhea
   cross-reference, or inferred from a database reaction with identical participants.
   Each suggestion states its evidence and links to ExPASy ENZYME and BRENDA.
2. **Enzyme candidates** — UniProtKB entries for the chosen EC, ranked with curated
   (SwissProt) entries first and organisms related to the production host promoted.
   This is the Selenzyme-style "reaction → clone this gene" step; see
   ``docs/enzyme_selection.md`` for why it is built on UniProt rather than the Selenzyme
   web service.
"""
from __future__ import annotations

from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QDesktopServices
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import (
    QComboBox, QDialog, QHBoxLayout, QHeaderView, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout)

from ...core import enzymes as ez
from .. import style
from ..widgets.busy import run_busy, was_cancelled


class EnzymeDialog(QDialog):
    def __init__(self, parent, reaction, suggestions: List[ez.ECSuggestion], *,
                 host_name: str = "", equation: str = ""):
        super().__init__(parent)
        self.setWindowTitle(f"Enzymes for {reaction.id}")
        self.resize(880, 600)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
        self._host_name = host_name
        self._suggestions = suggestions

        v = QVBoxLayout(self)
        head = QLabel(f"<b>{reaction.id}</b>"
                      + (f" — {reaction.name}" if reaction.name and
                         reaction.name != reaction.id else "")
                      + (f"<br><span style='color:#5f6368'>{equation}</span>"
                         if equation else ""))
        head.setWordWrap(True)
        head.setTextFormat(Qt.RichText)
        v.addWidget(head)

        if not suggestions:
            msg = QLabel(
                "No EC number could be found or inferred for this reaction.<br><br>"
                "<span style='color:#5f6368'>It carries no EC annotation, no KEGG/Rhea "
                "cross-reference to look one up from, and no reaction with identical "
                "participants in the loaded databases has one either. For rule-generated "
                "(RetroRules) chemistry this is expected — the transformation may not "
                "correspond to a classified enzyme at all.</span>")
            msg.setWordWrap(True)
            msg.setTextFormat(Qt.RichText)
            v.addWidget(msg)
        else:
            v.addWidget(QLabel("<b>Suggested EC numbers</b> (best evidence first):"))
            self.ec_table = QTableWidget(len(suggestions), 4)
            self.ec_table.setHorizontalHeaderLabels(
                ["EC number", "Evidence", "Source", "Links"])
            self.ec_table.verticalHeader().setVisible(False)
            self.ec_table.setEditTriggers(QTableWidget.NoEditTriggers)
            self.ec_table.setSelectionBehavior(QTableWidget.SelectRows)
            for r, s in enumerate(suggestions):
                self.ec_table.setItem(r, 0, QTableWidgetItem(s.ec))
                conf = {"annotated": "annotated ✓", "cross-reference": "cross-reference",
                        "inferred": "inferred — verify"}.get(s.confidence, s.confidence)
                self.ec_table.setItem(r, 1, QTableWidgetItem(conf))
                item = QTableWidgetItem(s.source)
                if s.detail:
                    item.setToolTip(s.detail)
                self.ec_table.setItem(r, 2, item)
                link = QLabel(f"<a href='{s.url}'>ExPASy</a> · "
                              f"<a href='{s.brenda_url}'>BRENDA</a>")
                link.setTextFormat(Qt.RichText)
                link.setOpenExternalLinks(True)
                self.ec_table.setCellWidget(r, 3, link)
            hh = self.ec_table.horizontalHeader()
            hh.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
            hh.setSectionResizeMode(2, QHeaderView.Stretch)
            hh.setSectionResizeMode(3, QHeaderView.ResizeToContents)
            self.ec_table.selectRow(0)
            v.addWidget(self.ec_table, 1)

            note = QLabel(
                "<span style='color:#8a6d00'>“inferred” means a different reaction with "
                "the same participants carries this EC — check the transformation really "
                "is the same before relying on it.</span>")
            note.setWordWrap(True)
            note.setTextFormat(Qt.RichText)
            v.addWidget(note)

            row = QHBoxLayout()
            row.addWidget(QLabel("Find enzymes for:"))
            self.ec_combo = QComboBox()
            self.ec_combo.addItems([s.ec for s in suggestions])
            row.addWidget(self.ec_combo)
            self.find_btn = QPushButton("Find enzyme candidates (UniProt)…")
            self.find_btn.setObjectName("primary")
            self.find_btn.setToolTip(
                "List UniProt entries for this EC, curated entries first and organisms "
                "related to your host promoted — i.e. candidate genes to clone.")
            self.find_btn.clicked.connect(self._find_enzymes)
            row.addWidget(self.find_btn)
            row.addStretch(1)
            v.addLayout(row)

            self.enz_table = QTableWidget(0, 5)
            self.enz_table.setHorizontalHeaderLabels(
                ["Accession", "Protein", "Organism", "Status", "Why ranked here"])
            self.enz_table.verticalHeader().setVisible(False)
            self.enz_table.setEditTriggers(QTableWidget.NoEditTriggers)
            self.enz_table.setSelectionBehavior(QTableWidget.SelectRows)
            self.enz_table.itemDoubleClicked.connect(self._open_uniprot)
            eh = self.enz_table.horizontalHeader()
            eh.setSectionResizeMode(1, QHeaderView.Stretch)
            eh.setSectionResizeMode(2, QHeaderView.Stretch)
            v.addWidget(self.enz_table, 1)
            self._candidates: List[ez.EnzymeCandidate] = []

        # ---- Selenzyme-style search by REACTION SIMILARITY ---------------------------
        # This is the branch that matters when the reaction has no EC at all (every
        # RetroRules step): it finds enzymes by comparing the chemistry itself.
        self._reaction = reaction
        from ...core import preferences, selenzyme as sz
        if preferences.selenzyme_enabled() and sz.is_installed():
            v.addWidget(QLabel(
                "<b>Or search by reaction similarity</b> — finds enzymes even with no "
                "EC number, by comparing this reaction to every characterised reaction "
                "in Rhea:"))
            sz_row = QHBoxLayout()
            self.sz_btn = QPushButton("Find enzymes by reaction similarity…")
            self.sz_btn.clicked.connect(self._find_by_similarity)
            sz_row.addWidget(self.sz_btn)
            self.sz_note = QLabel("")
            self.sz_note.setStyleSheet(f"color: {style.TEXT_MUTED};")
            sz_row.addWidget(self.sz_note, 1)
            v.addLayout(sz_row)

            self.sim_table = QTableWidget(0, 6)
            self.sim_table.setHorizontalHeaderLabels(
                ["Similarity", "Accession", "Protein", "Organism", "EC", "Closest Rhea"])
            self.sim_table.verticalHeader().setVisible(False)
            self.sim_table.setEditTriggers(QTableWidget.NoEditTriggers)
            self.sim_table.setSelectionBehavior(QTableWidget.SelectRows)
            self.sim_table.itemDoubleClicked.connect(self._open_sim_uniprot)
            sh = self.sim_table.horizontalHeader()
            sh.setSectionResizeMode(2, QHeaderView.Stretch)
            sh.setSectionResizeMode(3, QHeaderView.Stretch)
            v.addWidget(self.sim_table, 1)
            self._sim_hits = []

        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        bottom = QHBoxLayout()
        bottom.addStretch(1)
        bottom.addWidget(close)
        v.addLayout(bottom)

    def _find_by_similarity(self) -> None:
        from ...core import selenzyme as sz
        from PySide6.QtWidgets import QMessageBox
        smi = sz.reaction_smiles_for(self._reaction)
        if not smi:
            QMessageBox.information(
                self, "No structures available",
                "This reaction's metabolites carry no structures (SMILES/InChI), so its "
                "chemistry cannot be fingerprinted.<br><br>"
                "<span style='color:#5f6368'>Genome-scale database metabolites often "
                "carry only identifiers. Rule-based (RetroRules) steps always have "
                "structures, so similarity search works there.</span>")
            return
        ok, res = run_busy(
            self, "Comparing this reaction to every characterised reaction in Rhea…",
            lambda: sz.find_enzymes(smi, host_name=self._host_name, limit=25),
            title="Reaction-similarity enzyme search", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                QMessageBox.warning(self, "Search failed", str(res))
            return
        self._sim_hits = res
        self.sim_table.setRowCount(0)
        for h in res:
            r = self.sim_table.rowCount()
            self.sim_table.insertRow(r)
            self.sim_table.setItem(r, 0, QTableWidgetItem(f"{h.similarity:.3f}"))
            self.sim_table.setItem(r, 1, QTableWidgetItem(h.accession))
            self.sim_table.setItem(r, 2, QTableWidgetItem(h.protein))
            self.sim_table.setItem(r, 3, QTableWidgetItem(h.organism))
            self.sim_table.setItem(r, 4, QTableWidgetItem(h.ec))
            item = QTableWidgetItem(h.rhea_id)
            item.setToolTip(h.why)
            self.sim_table.setItem(r, 5, item)
        if res:
            best = res[0].similarity
            self.sz_note.setText(
                f"{len(res)} candidate(s); best reaction similarity {best:.2f}"
                + ("  — a close match" if best >= 0.6 else
                   "  — only a loose match, treat with caution" if best < 0.35 else ""))
            self.sim_table.selectRow(0)
        else:
            self.sz_note.setText("No sufficiently similar characterised reaction found.")

    def _open_sim_uniprot(self, item) -> None:
        row = item.row()
        if 0 <= row < len(self._sim_hits):
            QDesktopServices.openUrl(QUrl(self._sim_hits[row].url))

    def _find_enzymes(self) -> None:
        ec = self.ec_combo.currentText().strip()
        if not ec:
            return
        ok, res = run_busy(
            self, f"Searching UniProt for enzymes with EC {ec}…",
            lambda: ez.enzyme_candidates(ec, host_name=self._host_name, limit=40),
            title="Enzyme candidates", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                from PySide6.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Lookup failed", str(res))
            return
        self._candidates = res
        self.enz_table.setRowCount(0)
        for c in res:
            r = self.enz_table.rowCount()
            self.enz_table.insertRow(r)
            self.enz_table.setItem(r, 0, QTableWidgetItem(c.accession))
            self.enz_table.setItem(r, 1, QTableWidgetItem(c.protein))
            self.enz_table.setItem(r, 2, QTableWidgetItem(c.organism))
            self.enz_table.setItem(r, 3, QTableWidgetItem(c.reviewed))
            self.enz_table.setItem(r, 4, QTableWidgetItem(c.why))
        if res:
            self.enz_table.selectRow(0)

    def _open_uniprot(self, item) -> None:
        row = item.row()
        if 0 <= row < len(self._candidates):
            QDesktopServices.openUrl(QUrl(self._candidates[row].url))
