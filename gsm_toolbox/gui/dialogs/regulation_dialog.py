"""Tools ▸ Regulation — manage rule sets, choose which is active, and edit them.

Rule sets are files, handled like reaction databases: loaded once, kept in a local
library, and available in later sessions. Nothing is active until the user picks
something, and that default is deliberate — a rule set encodes one organism's physiology,
so applying one automatically would attach, say, cyanobacterial carbon regulation to
whatever model happened to be open and report the result as if it were biology.

The list therefore leads with **organism** and **confidence**. Those two columns answer
the questions that decide whether a rule set can be trusted for the model at hand, and
burying either in a detail pane would make it easy to skip.
"""
from __future__ import annotations

import os
from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QFileDialog, QHBoxLayout, QHeaderView, QLabel,
    QMenu, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QTextEdit,
    QVBoxLayout)

from ...core import preferences as prefs
from ...core import regulation as reg
from ...core import rule_library as lib
from .. import style

_CONFIDENCE_COLOUR = {
    reg.MEASURED: "#188038",
    reg.INFERRED: "#b06000",
    reg.ASSUMED: "#c5221f",
}


class RegulationDialog(QDialog):
    """The rule-set library: import, choose, inspect, edit."""

    def __init__(self, parent, model=None):
        super().__init__(parent)
        self.setWindowTitle("Regulatory rule sets")
        self.resize(1020, 620)
        from ..widgets.dialog_util import clamp_to_screen

        self._model = model
        self._entries: List[lib.RulesetInfo] = []

        v = QVBoxLayout(self)
        intro = QLabel(
            "Rule sets are files you load, not built-in defaults. Nothing is applied "
            "until you make a set active — regulation describes one organism, and "
            "applying the wrong organism's rules produces confident, wrong answers.")
        intro.setWordWrap(True)
        intro.setTextFormat(Qt.RichText)
        v.addWidget(intro)

        self.enable = QCheckBox("Apply the active rule set during simulation")
        self.enable.setChecked(bool(prefs.get(prefs.REGULATION_ENABLED)))
        self.enable.setToolTip(
            "When off, every simulation behaves exactly as it did before regulation "
            "existed.")
        v.addWidget(self.enable)

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            ["Active", "Rule set", "Organism", "Rules", "Weakest confidence", "Added"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.currentCellChanged.connect(lambda *_: self._show_detail())
        self.table.itemDoubleClicked.connect(lambda _: self._make_active())
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        for c in (0, 2, 3, 4, 5):
            header.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        v.addWidget(self.table, 1)

        self.detail = QTextEdit()
        self.detail.setReadOnly(True)
        self.detail.setMaximumHeight(170)
        self.detail.setPlaceholderText(
            "Select a rule set to see what it contains and how well it fits the loaded "
            "model.")
        v.addWidget(self.detail)

        row = QHBoxLayout()
        self.load_btn = QPushButton("Load rule file…")
        self.load_btn.setToolTip(
            "Import a rule-set JSON. It is validated, then kept in your library for "
            "later sessions.")
        self.load_btn.clicked.connect(self._load_file)
        row.addWidget(self.load_btn)

        self.example_btn = QPushButton("Load stored ▾")
        self.example_btn.setToolTip(
            "Rule-set files present in your regulation folder, plus any shipped as "
            "examples. Nothing is applied until you make it active.")
        self.example_btn.clicked.connect(self._show_examples)
        row.addWidget(self.example_btn)

        self.new_btn = QPushButton("New rule set…")
        self.new_btn.clicked.connect(self._new_ruleset)
        row.addWidget(self.new_btn)

        self.edit_btn = QPushButton("Edit…")
        self.edit_btn.clicked.connect(self._edit)
        row.addWidget(self.edit_btn)

        self.active_btn = QPushButton("Make active")
        self.active_btn.clicked.connect(self._make_active)
        row.addWidget(self.active_btn)

        self.remove_btn = QPushButton("Remove")
        self.remove_btn.clicked.connect(self._remove)
        row.addWidget(self.remove_btn)
        row.addStretch(1)
        v.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._save)
        buttons.rejected.connect(self.reject)
        v.addWidget(buttons)

        self._reload()
        clamp_to_screen(self)

    # -- library --------------------------------------------------------------------
    def _reload(self, select_path: str = "") -> None:
        self._entries = lib.scan_library()
        active = os.path.abspath(lib.active_path() or "")
        self.table.setRowCount(len(self._entries))
        chosen_row = -1
        for row, info in enumerate(self._entries):
            is_active = os.path.abspath(info.path) == active
            if is_active or (select_path and
                             os.path.abspath(info.path) == os.path.abspath(select_path)):
                chosen_row = row
            mark = QTableWidgetItem("●" if is_active else "")
            mark.setTextAlignment(Qt.AlignCenter)
            if is_active:
                mark.setToolTip("The rule set simulations will use.")
            self.table.setItem(row, 0, mark)

            name = QTableWidgetItem(info.name)
            if not info.valid:
                name.setForeground(QColor("#c5221f"))
                name.setToolTip("This file has problems and cannot be used.")
            self.table.setItem(row, 1, name)
            self.table.setItem(row, 2, QTableWidgetItem(info.organism or "—"))
            self.table.setItem(
                row, 3, QTableWidgetItem(f"{info.n_enabled} of {info.n_rules}"))
            confidence = info.weakest_confidence()
            citem = QTableWidgetItem(confidence)
            citem.setForeground(QColor(_CONFIDENCE_COLOUR.get(confidence, "#000")))
            self.table.setItem(row, 4, citem)
            self.table.setItem(row, 5, QTableWidgetItem(info.imported_text()))

        if not self._entries:
            self.detail.setPlainText(
                "Your library is empty.\n\n"
                "Use 'Load rule file…' to import one, 'Load example' for the rule sets "
                "shipped with the toolbox, or 'New rule set…' to write your own against "
                "the loaded model.\n\n"
                "Until a rule set is active, every simulation runs unregulated — exactly "
                "as it did before this feature existed.")
        elif chosen_row >= 0:
            self.table.selectRow(chosen_row)
        else:
            self.table.selectRow(0)
        self._sync_buttons()

    def _sync_buttons(self) -> None:
        info = self._selected()
        self.edit_btn.setEnabled(info is not None and info.valid)
        self.active_btn.setEnabled(info is not None and info.valid)
        self.remove_btn.setEnabled(info is not None)
        self.new_btn.setEnabled(self._model is not None)
        self.new_btn.setToolTip(
            "Write rules against the loaded model." if self._model is not None
            else "Load a model first — rules are written against its reactions.")

    def _selected(self) -> Optional[lib.RulesetInfo]:
        row = self.table.currentRow()
        if 0 <= row < len(self._entries):
            return self._entries[row]
        return None

    def _show_detail(self) -> None:
        info = self._selected()
        self._sync_buttons()
        if info is None:
            return
        lines = [f"<b>{info.name}</b>"]
        if info.organism:
            lines.append(f"Organism: <b>{info.organism}</b>")
        if info.description:
            lines.append(info.description)
        if info.source:
            lines.append(f"<i>{info.source}</i>")
        lines.append(f"{info.n_enabled} of {info.n_rules} rule(s) enabled.")

        if info.problems:
            lines.append("<span style='color:#c5221f'><b>This file cannot be used:</b>"
                         "<br>• " + "<br>• ".join(info.problems[:10]) + "</span>")
        elif self._model is not None:
            try:
                ruleset = reg.load(info.path)
                report = lib.fit(ruleset, self._model)
                colour = "#c5221f" if report.inapplicable else style.TEXT_MUTED
                lines.append(f"<span style='color:{colour}'>{report.summary()}</span>")

                # How many rules can actually change the answer, and why the rest
                # cannot. A rule that fires while changing nothing is the most
                # misleading state possible, so it is stated here rather than left to
                # be discovered by comparing two runs that came out the same.
                audit = lib.effectiveness(ruleset, self._model)
                warnings = audit.warnings()
                tone = "#b06000" if warnings else "#188038"
                block = f"<span style='color:{tone}'><b>{audit.summary()}</b>"
                if warnings:
                    block += "<br>• " + "<br>• ".join(w for w in warnings[:8])
                    if len(warnings) > 8:
                        block += f"<br>… and {len(warnings) - 8} more"
                lines.append(block + "</span>")
            except Exception as exc:  # noqa: BLE001
                lines.append(f"<span style='color:#c5221f'>Could not check against the "
                             f"loaded model: {exc}</span>")
        else:
            lines.append(f"<span style='color:{style.TEXT_MUTED}'>Load a model to check "
                         "whether these rules apply to it.</span>")
        self.detail.setHtml("<br><br>".join(lines))

    # -- actions --------------------------------------------------------------------
    def _load_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load regulatory rule set", "", "Rule sets (*.json);;All files (*)")
        if not path:
            return
        self._import(path)

    def _show_examples(self) -> None:
        """List every rule-set file on disk: the library folder first, then examples.

        The library folder is the one place a user can drop a file by hand, so it is
        listed by filename — a set copied in outside the app appears here without
        needing to be imported through a file dialog.
        """
        menu = QMenu(self)
        stored = lib.scan_library()
        if stored:
            header = menu.addAction(f"In {lib.library_dir()}")
            header.setEnabled(False)
            for info in stored:
                label = f"    {os.path.basename(info.path)}"
                if not info.valid:
                    label += "   (unusable)"
                action = menu.addAction(
                    label, lambda checked=False, p=info.path: self._activate(p))
                action.setEnabled(info.valid)
            menu.addSeparator()

        examples = [(label, path) for label, path in lib.examples()
                    if not any(os.path.basename(path) == os.path.basename(i.path)
                               for i in stored)]
        if examples:
            header = menu.addAction("Shipped examples (imported on click)")
            header.setEnabled(False)
            for label, path in examples:
                menu.addAction(f"    {label}",
                               lambda checked=False, p=path: self._import(p))
        if not stored and not examples:
            menu.addAction("No rule-set files found").setEnabled(False)

        self._example_menu = menu           # keep alive (see main_window._build_actions)
        menu.exec(self.example_btn.mapToGlobal(self.example_btn.rect().bottomLeft()))

    def _activate(self, path: str) -> None:
        """Make a file that is already in the library the active rule set."""
        info = lib.inspect(path)
        if not info.valid:
            QMessageBox.warning(
                self, "Cannot use this rule set",
                f"'{os.path.basename(path)}' has problems:\n\n• "
                + "\n• ".join(info.problems[:12]))
            return
        lib.set_active(path)
        self._reload(select_path=path)

    def _import(self, path: str) -> None:
        info = lib.import_file(path)
        if not info.valid:
            QMessageBox.warning(
                self, "Rule set not loaded",
                f"'{os.path.basename(path)}' was not imported:\n\n• "
                + "\n• ".join(info.problems[:12]))
            return
        self._reload(select_path=info.path)
        QMessageBox.information(
            self, "Rule set loaded",
            f"'{info.name}' added with {info.n_rules} rule(s).\n\n"
            "It is stored in your library and will be here next session. Use "
            "'Make active' to apply it.")

    def _new_ruleset(self) -> None:
        from .rule_editor_dialog import RuleEditorDialog

        empty = reg.RuleSet(name="New rule set")
        dialog = RuleEditorDialog(self, empty, self._model, "")
        if dialog.exec() and dialog.saved_path:
            self._reload(select_path=dialog.saved_path)

    def _edit(self) -> None:
        info = self._selected()
        if info is None:
            return
        from .rule_editor_dialog import RuleEditorDialog

        try:
            ruleset = reg.load(info.path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Cannot open", str(exc))
            return
        dialog = RuleEditorDialog(self, ruleset, self._model, info.path)
        if dialog.exec() and dialog.saved_path:
            self._reload(select_path=dialog.saved_path)

    def _make_active(self) -> None:
        info = self._selected()
        if info is None or not info.valid:
            return
        lib.set_active(info.path)
        self._reload(select_path=info.path)

    def _remove(self) -> None:
        info = self._selected()
        if info is None:
            return
        if QMessageBox.question(
                self, "Remove rule set",
                f"Remove '{info.name}' from your library?\n\nThe original file you "
                "imported is not affected.") != QMessageBox.Yes:
            return
        lib.remove(info.path)
        self._reload()

    def _save(self) -> None:
        prefs.set(prefs.REGULATION_ENABLED, bool(self.enable.isChecked()))
        self.accept()
