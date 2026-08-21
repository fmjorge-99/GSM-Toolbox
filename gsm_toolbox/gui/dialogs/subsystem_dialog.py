"""Review detected subsystems before they are written onto the model.

Detection is inference, not measurement. It reads reaction ids, EC numbers, pathway
intermediates and — as a last resort — reaction names, and any of those can be wrong for a
particular model. Writing the result straight onto the reactions would replace an honest
"unknown" with a confident label the user never saw, and every downstream feature that
groups by subsystem would then inherit the mistake silently.

So the flow is: detect → **review here** → apply. Nothing touches the model until the user
accepts.

The evidence behind each assignment is a column, not a tooltip. An assignment matched on
EC 5.3.1.9 and one matched because the reaction name contains "kinase" are not equally
trustworthy, and a reviewer who cannot see the difference has no basis for deciding which
rows to check.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QHBoxLayout,
    QHeaderView, QInputDialog, QLabel, QListWidget, QListWidgetItem, QMessageBox,
    QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from ...core import subsystems as sub
from .. import style

#: How strongly to flag each evidence level in the table.
_EVIDENCE_COLOUR = {
    sub.STRUCTURE: "#188038",
    sub.BY_ID: "#188038",
    sub.BY_EC: "#0b7285",
    sub.BY_METABOLITE: "#b06000",
    sub.BY_NAME: "#c5221f",
}

#: The pseudo-pathway holding reactions detection could not place.
_UNASSIGNED = "— Unassigned —"
#: The pseudo-pathway holding reactions that already carried a curated subsystem.
_EXISTING = "— Already annotated —"


class SubsystemDialog(QDialog):
    """Inspect, adjust and accept automatically detected subsystems."""

    def __init__(self, parent, model, report: sub.DetectionReport):
        super().__init__(parent)
        self.setWindowTitle("Detected subsystems")
        self.resize(1100, 660)
        from ..widgets.dialog_util import clamp_to_screen

        self._model = model
        self._report = report
        #: reaction id → pathway name. The single source of truth the dialog edits.
        self._assign: Dict[str, str] = {
            rid: a.pathway for rid, a in report.assignments.items()}
        self._evidence: Dict[str, str] = {
            rid: a.evidence_text() for rid, a in report.assignments.items()}
        self._evidence_kind: Dict[str, str] = {
            rid: a.evidence for rid, a in report.assignments.items()}
        #: pathway → apply it? Detection proposes; the user decides which proposals to
        #: take. All start ticked, so accepting everything stays one click.
        self._selected: Dict[str, bool] = {
            a.pathway: True for a in report.assignments.values()}
        self._loading = False

        outer = QVBoxLayout(self)
        self.summary = QLabel(report.summary())
        self.summary.setWordWrap(True)
        outer.addWidget(self.summary)

        note = QLabel(
            "These assignments are <b>inferred</b>, not read from the model. Review them "
            "— especially rows marked <span style='color:#c5221f'>reaction name</span>, "
            "which is the weakest evidence — then accept to write them onto the "
            "reactions.")
        note.setWordWrap(True)
        note.setTextFormat(Qt.RichText)
        outer.addWidget(note)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._pathway_panel())
        splitter.addWidget(self._reaction_panel())
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([330, 770])
        outer.addWidget(splitter, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.apply_btn = buttons.addButton("Apply to model",
                                           QDialogButtonBox.AcceptRole)
        self.apply_btn.clicked.connect(self._accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

        self._reload_pathways()
        clamp_to_screen(self)

    # -- panels ---------------------------------------------------------------------
    def _pathway_panel(self) -> QWidget:
        panel = QWidget()
        v = QVBoxLayout(panel)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(QLabel("<b>Detected pathways</b> — tick the ones to apply"))
        self.pathways = QListWidget()
        self.pathways.currentRowChanged.connect(lambda _: self._show_reactions())
        self.pathways.itemChanged.connect(self._on_pathway_toggled)
        v.addWidget(self.pathways, 1)

        select_row = QHBoxLayout()
        for label, state in (("Select all", True), ("Select none", False)):
            button = QPushButton(label)
            button.clicked.connect(lambda _=False, s=state: self._set_all_checked(s))
            select_row.addWidget(button)
        v.addLayout(select_row)

        rename = QPushButton("Rename pathway…")
        rename.setToolTip(
            "Use the naming convention your other models follow — the name written onto "
            "the reactions is whatever appears here.")
        rename.clicked.connect(self._rename_pathway)
        v.addWidget(rename)

        self.detail = QLabel("")
        self.detail.setWordWrap(True)
        self.detail.setStyleSheet(f"color:{style.TEXT_MUTED};")
        v.addWidget(self.detail)
        return panel

    def _reaction_panel(self) -> QWidget:
        panel = QWidget()
        v = QVBoxLayout(panel)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(QLabel("<b>Reactions</b>"))

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["Reaction", "Name", "Why it was assigned", "Subsystem to write"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        v.addWidget(self.table, 1)

        row = QHBoxLayout()
        add = QPushButton("Add reactions…")
        add.setToolTip("Assign more reactions to this pathway, from anywhere in the "
                       "model.")
        add.clicked.connect(self._add_reactions)
        row.addWidget(add)

        move = QPushButton("Move to…")
        move.setToolTip("Reassign the selected reactions to another pathway.")
        move.clicked.connect(self._move_selected)
        row.addWidget(move)

        remove = QPushButton("Remove from pathway")
        remove.setToolTip("Leave the selected reactions unassigned.")
        remove.clicked.connect(self._remove_selected)
        row.addWidget(remove)
        row.addStretch(1)
        v.addLayout(row)
        return panel

    # -- pathway list ---------------------------------------------------------------
    def _pathway_names(self) -> List[str]:
        """Pathways currently holding at least one reaction, in catalogue order."""
        used = set(self._assign.values())
        ordered = [p.name for p in sub.CATALOGUE if p.name in used]
        extra = sorted(used - set(ordered))       # user-renamed or added
        return ordered + extra

    def _reactions_in(self, pathway: str) -> List[str]:
        if pathway == _UNASSIGNED:
            assigned = set(self._assign)
            return [r for r in self._report.unassigned if r not in assigned]
        if pathway == _EXISTING:
            return sorted(self._report.preserved)
        return sorted(rid for rid, name in self._assign.items() if name == pathway)

    def _reload_pathways(self, keep: str = "") -> None:
        keep = keep or self._current_pathway() or ""
        self._loading = True
        self.pathways.clear()
        for name in self._pathway_names():
            count = len(self._reactions_in(name))
            item = QListWidgetItem(f"{name}  ({count})")
            item.setData(Qt.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if self._selected.get(name, True)
                               else Qt.Unchecked)
            item.setToolTip("Only ticked pathways are written to the model.")
            self.pathways.addItem(item)

        for pseudo, tip in ((_UNASSIGNED, "Detection could not place these."),
                            (_EXISTING, "Already annotated in the model; untouched.")):
            entries = self._reactions_in(pseudo)
            if not entries:
                continue
            item = QListWidgetItem(f"{pseudo}  ({len(entries)})")
            item.setData(Qt.UserRole, pseudo)
            item.setForeground(QColor(style.TEXT_MUTED))
            item.setToolTip(tip)
            self.pathways.addItem(item)
        self._loading = False

        target = 0
        for row in range(self.pathways.count()):
            if self.pathways.item(row).data(Qt.UserRole) == keep:
                target = row
                break
        self.pathways.setCurrentRow(target)
        self._update_summary()

    def _current_pathway(self) -> Optional[str]:
        item = self.pathways.currentItem()
        return item.data(Qt.UserRole) if item else None

    def _on_pathway_toggled(self, item: QListWidgetItem) -> None:
        if self._loading:
            return
        name = item.data(Qt.UserRole)
        if name in (_UNASSIGNED, _EXISTING):
            return
        self._selected[name] = item.checkState() == Qt.Checked
        self._update_summary()

    def _set_all_checked(self, state: bool) -> None:
        for name in self._pathway_names():
            self._selected[name] = state
        self._reload_pathways()

    def _update_summary(self) -> None:
        chosen = self.assignments()
        pathways = {p for p in chosen.values()}
        skipped = len(self._assign) - len(chosen)
        weak = sum(1 for rid in chosen
                   if self._evidence_kind.get(rid) == sub.BY_NAME)
        text = (f"{len(chosen)} reaction(s) will be annotated across "
                f"{len(pathways)} subsystem(s).")
        if skipped:
            text += f" {skipped} in unticked pathways will be skipped."
        if self._report.preserved:
            text += (f" {len(self._report.preserved)} already annotated and left "
                     f"unchanged.")
        remaining = len(self._reactions_in(_UNASSIGNED))
        if remaining:
            text += f" {remaining} could not be placed."
        if weak:
            text += f" {weak} rest on a name match only."
        self.summary.setText(text)
        self.apply_btn.setEnabled(bool(chosen))

    def _rename_pathway(self) -> None:
        current = self._current_pathway()
        if not current or current in (_UNASSIGNED, _EXISTING):
            return
        name, ok = QInputDialog.getText(self, "Rename pathway", "Subsystem name:",
                                        text=current)
        name = (name or "").strip()
        if not ok or not name or name == current:
            return
        for rid, pathway in list(self._assign.items()):
            if pathway == current:
                self._assign[rid] = name
        self._reload_pathways(keep=name)

    # -- reaction table -------------------------------------------------------------
    def _show_reactions(self) -> None:
        if self._loading:
            return
        pathway = self._current_pathway()
        if pathway is None:
            return
        reactions = self._reactions_in(pathway)
        self.table.setRowCount(len(reactions))
        for row, rid in enumerate(reactions):
            name = ""
            if self._model.reactions.has_id(rid):
                name = self._model.reactions.get_by_id(rid).name or ""
            self.table.setItem(row, 0, QTableWidgetItem(rid))
            self.table.setItem(row, 1, QTableWidgetItem(name))

            if pathway == _EXISTING:
                evidence = "declared by the model"
                kind = sub.STRUCTURE
                written = self._report.preserved.get(rid, "")
            elif pathway == _UNASSIGNED:
                evidence = "no signal matched"
                kind = sub.BY_NAME
                written = "(none)"
            else:
                evidence = self._evidence.get(rid, "chosen by hand")
                kind = self._evidence_kind.get(rid, sub.STRUCTURE)
                written = pathway
            item = QTableWidgetItem(evidence)
            item.setForeground(QColor(_EVIDENCE_COLOUR.get(kind, "#000")))
            self.table.setItem(row, 2, item)
            self.table.setItem(row, 3, QTableWidgetItem(written))

        described = next((p for p in sub.CATALOGUE if p.name == pathway), None)
        if described:
            self.detail.setText(described.description)
        elif pathway == _UNASSIGNED:
            self.detail.setText(
                "Nothing matched these. Specialised pathways outside central metabolism "
                "are expected here — leaving them blank is more useful than a wrong "
                "label.")
        elif pathway == _EXISTING:
            self.detail.setText(
                "These already carry a subsystem. Curated annotation beats inference, so "
                "they are left alone.")
        else:
            self.detail.setText("")

    def _selected_reactions(self) -> List[str]:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        return [self.table.item(r, 0).text() for r in rows]

    def _add_reactions(self) -> None:
        pathway = self._current_pathway()
        if not pathway or pathway == _EXISTING:
            QMessageBox.information(
                self, "Pick a pathway",
                "Select the pathway to add reactions to first.")
            return
        from .reaction_browser_dialog import ReactionBrowserDialog

        picked = ReactionBrowserDialog.pick(
            self, self._model, self._reactions_in(pathway),
            title=f"Reactions in '{pathway}'",
            prompt="Everything in the model is listed. Filter by subsystem or search, "
                   "and use <b>Add all shown</b> to take a whole group at once.")
        if picked is None:
            return
        chosen = set(picked)
        # Anything dropped from the selection leaves the pathway; anything added joins it,
        # moving out of whatever pathway held it before.
        for rid in self._reactions_in(pathway):
            if rid not in chosen:
                self._assign.pop(rid, None)
        for rid in chosen:
            if self._assign.get(rid) != pathway:
                self._assign[rid] = pathway
                self._evidence[rid] = "chosen by hand"
                self._evidence_kind[rid] = sub.STRUCTURE
        self._reload_pathways(keep=pathway)
        self._show_reactions()

    def _move_selected(self) -> None:
        reactions = self._selected_reactions()
        if not reactions:
            return
        options = self._pathway_names() + [p.name for p in sub.CATALOGUE
                                           if p.name not in self._pathway_names()]
        name, ok = QInputDialog.getItem(self, "Move reactions",
                                        f"Move {len(reactions)} reaction(s) to:",
                                        options, 0, True)
        if not ok or not name.strip():
            return
        for rid in reactions:
            self._assign[rid] = name.strip()
            self._evidence[rid] = "chosen by hand"
            self._evidence_kind[rid] = sub.STRUCTURE
        self._reload_pathways(keep=name.strip())
        self._show_reactions()

    def _remove_selected(self) -> None:
        pathway = self._current_pathway()
        if pathway == _EXISTING:
            QMessageBox.information(
                self, "Not editable here",
                "These reactions already carry a subsystem from the model and are not "
                "being changed. Re-run detection with 'replace existing' to revisit "
                "them.")
            return
        for rid in self._selected_reactions():
            self._assign.pop(rid, None)
        self._reload_pathways(keep=pathway or "")
        self._show_reactions()

    # -- result ---------------------------------------------------------------------
    def assignments(self) -> Dict[str, str]:
        """Only the reactions in ticked pathways — the rest are left untouched."""
        return {rid: pathway for rid, pathway in self._assign.items()
                if self._selected.get(pathway, True)}

    def _accept(self) -> None:
        if not self.assignments():
            QMessageBox.information(
                self, "Nothing selected",
                "Tick at least one pathway to write onto the model.")
            return
        self.accept()


class SubsystemOptionsDialog(QDialog):
    """Asked before detection runs — the two choices that change what it does."""

    def __init__(self, parent, model):
        super().__init__(parent)
        self.setWindowTitle("Detect subsystems")
        self.resize(560, 300)
        from ..widgets.dialog_util import clamp_to_screen

        existing = sub.existing_subsystems(model)
        v = QVBoxLayout(self)

        intro = QLabel(
            "Assigns each reaction to a central pathway, using reaction ids, EC "
            "numbers, pathway intermediates and stoichiometry.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        if existing:
            state = QLabel(
                f"<b>This model already declares {len(existing)} subsystem(s)</b> across "
                f"{sum(existing.values())} reactions.")
        else:
            state = QLabel("<b>This model declares no subsystems.</b>")
        state.setWordWrap(True)
        state.setTextFormat(Qt.RichText)
        v.addWidget(state)

        self.overwrite = QCheckBox("Replace subsystems the model already declares")
        self.overwrite.setChecked(False)
        self.overwrite.setEnabled(bool(existing))
        self.overwrite.setToolTip(
            "Off by default: a curated annotation is better evidence than anything "
            "inferred here, so existing ones are left alone.")
        v.addWidget(self.overwrite)

        v.addWidget(QLabel("Weakest evidence to accept:"))
        self.evidence = QComboBox()
        self.evidence.addItem("Reaction name — most coverage, needs review",
                              sub.BY_NAME)
        self.evidence.addItem("Pathway intermediates and stronger", sub.BY_METABOLITE)
        self.evidence.addItem("EC numbers and stronger", sub.BY_EC)
        self.evidence.addItem("Reaction ids and stoichiometry only — most conservative",
                              sub.BY_ID)
        self.evidence.setToolTip(
            "Detection layers five signals. This sets how weak a signal may be before "
            "an assignment is dropped; everything accepted is shown with its evidence "
            "for review either way.")
        v.addWidget(self.evidence)
        v.addStretch(1)

        buttons = QDialogButtonBox(QDialogButtonBox.Cancel)
        buttons.addButton("Detect", QDialogButtonBox.AcceptRole
                          ).clicked.connect(self.accept)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)
        clamp_to_screen(self)

    def values(self) -> dict:
        return {"overwrite": self.overwrite.isChecked(),
                "minimum_evidence": self.evidence.currentData()}
