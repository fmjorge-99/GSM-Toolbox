"""Explore heterologous chemistry upstream of a designed route's entry point.

A design ends where it meets native metabolism, and that hand-over is where the yield is
usually decided — yet it is the part a route report says least about. This dialog opens it
up: for each compound the route draws from the host it lists the database reactions that
could produce that compound, marked by whether the host can actually run them today.

Two situations bring a user here:

* the entry compound is **idle** — present in the model but carrying no flux — so the
  route cannot run until something supplies it (the lactaldehyde case);
* the entry compound does carry flux, but a heterologous step feeding it would give the
  whole route more carbon.

Selected reactions can be added to the route, so the design grows upstream in place.
"""
from __future__ import annotations

from typing import List

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QDialog, QHBoxLayout, QHeaderView, QLabel,
    QPushButton, QTableWidget, QTableWidgetItem, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout)

from .. import style

_READY_COLOUR = "#188038"
_BLOCKED_COLOUR = "#b06000"
_MISSING_COLOUR = "#c5221f"


class UpstreamDialog(QDialog):
    """Lists upstream candidates per entry metabolite; emits the reactions to add."""

    reactions_chosen = Signal(list)

    def __init__(self, parent, report):
        super().__init__(parent)
        self.setWindowTitle("Upstream of this route")
        self.resize(940, 620)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
        self._report = report

        v = QVBoxLayout(self)
        head = QLabel(report.headline())
        head.setWordWrap(True)
        head.setTextFormat(Qt.RichText)
        head.setStyleSheet("padding:8px; background:#F1F3F4; border-radius:4px;")
        v.addWidget(head)

        v.addWidget(_muted(
            "Each entry compound is listed with the reactions that could make it. "
            "<b>Green</b> means every substrate is already produced by your host, so the "
            "step would run as soon as the enzyme is expressed. <b>Amber</b> means the "
            "substrates exist but are themselves idle — that step needs its own supply. "
            "<b>Red</b> means a substrate is absent from the host altogether."))

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(
            ["Reaction", "Status", "Equation", "EC"])
        self.tree.setSelectionMode(QAbstractItemView.NoSelection)
        self.tree.setAlternatingRowColors(True)
        hh = self.tree.header()
        hh.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.Stretch)
        hh.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        v.addWidget(self.tree, 1)
        self._populate()

        self.ready_only = QCheckBox("Show only steps the host could run today")
        self.ready_only.setToolTip(
            "Hide candidates whose own substrates are missing or idle. Useful for "
            "picking a step that works without a second round of engineering.")
        self.ready_only.toggled.connect(self._apply_filter)
        v.addWidget(self.ready_only)

        v.addWidget(_muted(report.context_summary))

        row = QHBoxLayout()
        self.add_btn = QPushButton("Add ticked reactions to the route")
        self.add_btn.setObjectName("primary")
        self.add_btn.setToolTip(
            "Extend the current pathway with the reactions you ticked, then re-analyse.")
        self.add_btn.clicked.connect(self._emit_choice)
        row.addWidget(self.add_btn)
        row.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        row.addWidget(close)
        v.addLayout(row)

    # -- contents ------------------------------------------------------------------
    def _populate(self) -> None:
        self.tree.clear()
        for mid, name in self._report.entry_metabolites.items():
            idle = mid in self._report.idle_entries
            label = f"{name}  ({mid})"
            if idle:
                label += "   — carries no flux"
            parent = QTreeWidgetItem([label, "entry point", "", ""])
            parent.setFirstColumnSpanned(False)
            if idle:
                parent.setForeground(0, _colour(_BLOCKED_COLOUR))
            parent.setToolTip(
                0, "The compound where this route draws on native metabolism.")
            self.tree.addTopLevelItem(parent)
            for cand in self._report.candidates.get(mid, []):
                if cand.missing_substrates:
                    colour, status = _MISSING_COLOUR, "substrate missing"
                elif cand.idle_substrates:
                    colour, status = _BLOCKED_COLOUR, "substrate idle"
                else:
                    colour, status = _READY_COLOUR, "ready"
                child = QTreeWidgetItem([
                    cand.reaction_id, status, cand.equation,
                    ", ".join(cand.ec_numbers)])
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
                child.setCheckState(0, Qt.Unchecked)
                child.setForeground(1, _colour(colour))
                tip = cand.verdict()
                if not cand.balance_checkable:
                    tip += " · balance could not be checked"
                elif not cand.balanced:
                    tip += " · mass/charge UNBALANCED"
                child.setToolTip(0, tip)
                child.setToolTip(1, tip)
                child.setData(0, Qt.UserRole, cand.reaction_id)
                child.setData(1, Qt.UserRole, bool(cand.ready))
                parent.addChild(child)
            parent.setExpanded(True)

    def _apply_filter(self, ready_only: bool) -> None:
        for i in range(self.tree.topLevelItemCount()):
            parent = self.tree.topLevelItem(i)
            for j in range(parent.childCount()):
                child = parent.child(j)
                hide = ready_only and not bool(child.data(1, Qt.UserRole))
                child.setHidden(hide)

    def _emit_choice(self) -> None:
        chosen: List[str] = []
        for i in range(self.tree.topLevelItemCount()):
            parent = self.tree.topLevelItem(i)
            for j in range(parent.childCount()):
                child = parent.child(j)
                if child.checkState(0) == Qt.Checked:
                    chosen.append(str(child.data(0, Qt.UserRole)))
        if chosen:
            self.reactions_chosen.emit(chosen)
            self.accept()


def _muted(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setTextFormat(Qt.RichText)
    lbl.setStyleSheet(f"color:{style.TEXT_MUTED};")
    return lbl


def _colour(hex_colour: str):
    from PySide6.QtGui import QColor
    return QColor(hex_colour)
