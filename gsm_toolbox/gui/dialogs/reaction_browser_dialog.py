"""Pick reactions from a model: everything available on the left, the chosen set right.

Typing reaction ids into a text field only works if you already know them. A genome-scale
model has thousands, named in a convention the user may not share, so the field was in
practice unusable — this replaces it.

Used in three places, which is why it is a dialog rather than a widget buried in one panel:

* choosing which fluxes a condition scan reports,
* choosing which fluxes a time course records,
* choosing which reactions a regulatory rule acts on.

Selecting by **subsystem** is offered alongside individual reactions. A regulatory rule
usually targets a pathway ("the carbon-concentrating mechanism"), not a reaction picked
one at a time, and requiring the latter would make writing a rule far more tedious than
the biology it describes.
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QHBoxLayout, QLabel,
    QLineEdit, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout, QWidget)

from .. import style

#: Shown in the subsystem filter when a reaction declares none.
_NO_SUBSYSTEM = "(no subsystem)"


class ReactionBrowserDialog(QDialog):
    """Two-panel reaction chooser with search, subsystem filter and bulk moves."""

    def __init__(self, parent, model, selected: Sequence[str] = (), *,
                 title: str = "Select reactions",
                 prompt: str = "", allow_metabolites: bool = False,
                 restrict_to: Optional[Iterable[str]] = None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 560)
        from ..widgets.dialog_util import clamp_to_screen

        self._model = model
        self._allow_metabolites = allow_metabolites
        self._restrict = set(restrict_to) if restrict_to is not None else None
        self._chosen: List[str] = [r for r in selected if r]

        outer = QVBoxLayout(self)
        if prompt:
            hint = QLabel(prompt)
            hint.setWordWrap(True)
            hint.setTextFormat(Qt.RichText)
            outer.addWidget(hint)

        filters = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search id, name or gene…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._refill)
        filters.addWidget(self.search, 2)

        self.subsystem = QComboBox()
        self.subsystem.addItem("All subsystems", "")
        for name in self._subsystems():
            self.subsystem.addItem(name, name)
        self.subsystem.currentIndexChanged.connect(self._refill)
        filters.addWidget(self.subsystem, 1)
        outer.addLayout(filters)

        panels = QHBoxLayout()
        panels.addWidget(self._panel("Available", "available"), 1)

        middle = QVBoxLayout()
        middle.addStretch(1)
        self.add_btn = QPushButton("Add  ▶")
        self.add_btn.clicked.connect(self._add_selected)
        self.add_all_btn = QPushButton("Add all shown  ▶▶")
        self.add_all_btn.setToolTip(
            "Adds every reaction matching the current search and subsystem filter — the "
            "quick way to target a whole pathway.")
        self.add_all_btn.clicked.connect(self._add_all_shown)
        self.remove_btn = QPushButton("◀  Remove")
        self.remove_btn.clicked.connect(self._remove_selected)
        self.clear_btn = QPushButton("◀◀ Remove all")
        self.clear_btn.clicked.connect(self._remove_all)
        for b in (self.add_btn, self.add_all_btn, self.remove_btn, self.clear_btn):
            b.setMinimumWidth(130)
            middle.addWidget(b)
        middle.addStretch(1)
        panels.addLayout(middle)

        panels.addWidget(self._panel("Selected", "selected"), 1)
        outer.addLayout(panels, 1)

        self.count = QLabel("")
        self.count.setStyleSheet(f"color:{style.TEXT_MUTED};")
        outer.addWidget(self.count)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

        self._refill()
        self._refresh_selected()
        clamp_to_screen(self)

    # -- construction helpers -------------------------------------------------------
    def _panel(self, title: str, which: str) -> QWidget:
        box = QWidget()
        v = QVBoxLayout(box)
        v.setContentsMargins(0, 0, 0, 0)
        label = QLabel(f"<b>{title}</b>")
        label.setTextFormat(Qt.RichText)
        v.addWidget(label)
        listing = QListWidget()
        listing.setSelectionMode(QAbstractItemView.ExtendedSelection)
        listing.setAlternatingRowColors(True)
        if which == "available":
            listing.itemDoubleClicked.connect(lambda _: self._add_selected())
            self.available = listing
        else:
            listing.itemDoubleClicked.connect(lambda _: self._remove_selected())
            self.selected = listing
        v.addWidget(listing, 1)
        return box

    def _subsystems(self) -> List[str]:
        if self._allow_metabolites:
            return []
        names = {(r.subsystem or _NO_SUBSYSTEM).strip() or _NO_SUBSYSTEM
                 for r in self._model.reactions}
        return sorted(names)

    # -- data -----------------------------------------------------------------------
    def _entries(self) -> Iterable[tuple]:
        """``(id, display)`` for everything selectable."""
        from ...core.network_graph import clean_label

        if self._allow_metabolites:
            for met in self._model.metabolites:
                if self._restrict is not None and met.id not in self._restrict:
                    continue
                yield met.id, clean_label(f"{met.id} — {met.name or ''}")
        else:
            for rxn in self._model.reactions:
                if self._restrict is not None and rxn.id not in self._restrict:
                    continue
                subsystem = (rxn.subsystem or "").strip()
                suffix = f"  [{subsystem}]" if subsystem else ""
                yield rxn.id, clean_label(f"{rxn.id} — {rxn.name or ''}{suffix}")

    def _matches(self, rid: str, display: str) -> bool:
        text = self.search.text().strip().lower()
        if text and text not in display.lower():
            return False
        wanted = self.subsystem.currentData()
        if wanted and not self._allow_metabolites:
            if not self._model.reactions.has_id(rid):
                return False
            actual = (self._model.reactions.get_by_id(rid).subsystem
                      or _NO_SUBSYSTEM).strip() or _NO_SUBSYSTEM
            if actual != wanted:
                return False
        return True

    def _refill(self) -> None:
        chosen = set(self._chosen)
        self.available.clear()
        shown = 0
        for rid, display in self._entries():
            if rid in chosen or not self._matches(rid, display):
                continue
            item = QListWidgetItem(display)
            item.setData(Qt.UserRole, rid)
            self.available.addItem(item)
            shown += 1
            if shown >= 3000:      # a full universal model would freeze the list
                item = QListWidgetItem("… refine the search to see more")
                item.setFlags(Qt.NoItemFlags)
                self.available.addItem(item)
                break
        self._update_count()

    def _refresh_selected(self) -> None:
        self.selected.clear()
        from ...core.network_graph import clean_label

        for rid in self._chosen:
            if not self._allow_metabolites and self._model.reactions.has_id(rid):
                rxn = self._model.reactions.get_by_id(rid)
                display = clean_label(f"{rid} — {rxn.name or ''}")
            elif self._allow_metabolites and self._model.metabolites.has_id(rid):
                met = self._model.metabolites.get_by_id(rid)
                display = clean_label(f"{rid} — {met.name or ''}")
            else:
                display = f"{clean_label(rid)}   (not in this model)"
            item = QListWidgetItem(display)
            item.setData(Qt.UserRole, rid)
            if "not in this model" in display:
                from PySide6.QtGui import QColor
                item.setForeground(QColor("#c5221f"))
                item.setToolTip("Kept, but it will have no effect on the loaded model.")
            self.selected.addItem(item)
        self._update_count()

    def _update_count(self) -> None:
        self.count.setText(f"{len(self._chosen)} selected · "
                           f"{self.available.count()} shown of "
                           f"{sum(1 for _ in self._entries())} in the model")

    # -- moves ----------------------------------------------------------------------
    def _add_selected(self) -> None:
        for item in self.available.selectedItems():
            rid = item.data(Qt.UserRole)
            if rid and rid not in self._chosen:
                self._chosen.append(rid)
        self._refill()
        self._refresh_selected()

    def _add_all_shown(self) -> None:
        for row in range(self.available.count()):
            rid = self.available.item(row).data(Qt.UserRole)
            if rid and rid not in self._chosen:
                self._chosen.append(rid)
        self._refill()
        self._refresh_selected()

    def _remove_selected(self) -> None:
        drop = {item.data(Qt.UserRole) for item in self.selected.selectedItems()}
        self._chosen = [r for r in self._chosen if r not in drop]
        self._refill()
        self._refresh_selected()

    def _remove_all(self) -> None:
        self._chosen = []
        self._refill()
        self._refresh_selected()

    # -- result ---------------------------------------------------------------------
    def chosen(self) -> List[str]:
        return list(self._chosen)

    @staticmethod
    def pick(parent, model, selected: Sequence[str] = (), **kwargs) -> Optional[List[str]]:
        """Show the dialog; returns the chosen ids, or None if cancelled."""
        dialog = ReactionBrowserDialog(parent, model, selected, **kwargs)
        if dialog.exec():
            return dialog.chosen()
        return None
