"""Side-by-side comparison of the routes currently open as result tabs (VI.9).

Each route lives in its own tab, so choosing between four alternatives means clicking
through four tabs and remembering four numbers. This puts them in one table, with the
figures that actually decide the choice — most importantly **carbon yield**, which is
the only flux-derived number that is comparable across routes and across targets: a
raw production flux of 12 mmol/gDW/h means nothing next to one of 3 until you know how
much carbon each consumed to get there.
"""
from __future__ import annotations

from typing import List

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QHeaderView, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout)

from .. import style

_HEADERS = ["Route", "Steps", "Carbon yield", "Production flux", "Chemistry",
            "Verdict"]


class CompareRoutesDialog(QDialog):
    """Read-only table of every open route, sorted by carbon yield."""

    def __init__(self, parent, results: List, titles: List[str]):
        super().__init__(parent)
        self.setWindowTitle("Compare routes")
        self.resize(900, 460)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

        v = QVBoxLayout(self)
        v.addWidget(QLabel(
            "<b>Rank on carbon yield</b> — mol carbon in the product per mol carbon "
            "consumed. Unlike raw flux, it is comparable between routes and targets."))

        rows = self._collect(results, titles)
        self.table = QTableWidget(len(rows), len(_HEADERS))
        self.table.setHorizontalHeaderLabels(_HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        for r, row in enumerate(rows):
            for c, (text, tip, colour, align) in enumerate(row["cells"]):
                item = QTableWidgetItem(text)
                if tip:
                    item.setToolTip(tip)
                if colour:
                    item.setForeground(QColor(colour))
                item.setTextAlignment(align)
                self.table.setItem(r, c, item)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        for c in range(1, len(_HEADERS) - 1):
            hh.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(len(_HEADERS) - 1, QHeaderView.Stretch)
        v.addWidget(self.table, 1)

        note = QLabel(
            "“—” in Carbon yield means the figure could not be computed: the route "
            "carries no flux, or its flux is rule-derived and therefore indicative "
            "only. Those routes are listed last and cannot be ranked on yield.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color:{style.TEXT_MUTED};")
        v.addWidget(note)

        btns = QHBoxLayout()
        btns.addStretch(1)
        close = QPushButton("Close")
        close.setObjectName("primary")
        close.clicked.connect(self.accept)
        btns.addWidget(close)
        v.addLayout(btns)

    # -- assembly -------------------------------------------------------------------
    def _collect(self, results: List, titles: List[str]) -> List[dict]:
        rows = []
        for result, title in zip(results, titles):
            cy = _f(getattr(result, "carbon_yield", float("nan")))
            flux = _f(result.production_flux)
            indicative = bool(getattr(result, "flux_is_indicative", False))
            rows.append({
                "sort": (-cy if cy == cy and not indicative else 1.0, title),
                "cells": [
                    (title, result.target, "", Qt.AlignLeft | Qt.AlignVCenter),
                    (str(len(result.reaction_ids)), "Heterologous steps to add", "",
                     Qt.AlignCenter),
                    _yield_cell(cy, indicative,
                                getattr(result, "carbon_yield_note", "") or ""),
                    _flux_cell(flux, indicative),
                    _chemistry_cell(result),
                    _verdict_cell(result),
                ],
            })
        rows.sort(key=lambda r: r["sort"])
        return rows


def _f(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _yield_cell(cy: float, indicative: bool, note: str = ""):
    align = Qt.AlignRight | Qt.AlignVCenter
    if indicative:
        return ("—", "Rule-derived route: its flux is indicative, so no meaningful "
                "carbon yield can be derived from it.", style.TEXT_MUTED, align)
    if cy != cy:
        # Say which of the several possible causes it actually was — "no flux" and
        # "the product has no formula" call for completely different next steps.
        tip = (f"Carbon yield not computable — {note}."
               if note else
               "No carbon yield — the route carries no flux, or the carbon "
               "content of the substrate or product could not be determined.")
        return ("—", tip, style.TEXT_MUTED, align)
    # Anything above roughly half the theoretical carbon is a strong route; below a
    # tenth the product is a trace by-product of whatever else the cell is doing.
    colour = "#188038" if cy >= 0.30 else ("#b06000" if cy >= 0.10 else "#c5221f")
    return (f"{cy * 100:.1f}%", "mol C in product / mol C consumed", colour, align)


def _flux_cell(flux: float, indicative: bool):
    align = Qt.AlignRight | Qt.AlignVCenter
    if flux != flux:
        return ("—", "", style.TEXT_MUTED, align)
    if indicative:
        return ("carries flux" if flux > 1e-9 else "no flux",
                "Rule-based route: only the yes/no is meaningful.",
                style.TEXT_MUTED, align)
    return (f"{flux:.4g}", "mmol gDW⁻¹ h⁻¹", "", align)


def _chemistry_cell(result):
    align = Qt.AlignLeft | Qt.AlignVCenter
    bits, tips = [], []
    if getattr(result, "isomer_warnings", None):
        bits.append("isomer ✗")
        tips.extend(result.isomer_warnings[:3])
    elif getattr(result, "isomer_checked", False):
        bits.append("isomer ✓")
    unverified = list(getattr(result, "unverified_steps", None) or [])
    n_steps = len(getattr(result, "reaction_ids", ()) or ())
    if not getattr(result, "balanced", True):
        bits.append("balance ✗")
        tips.append("At least one step is mass/charge-unbalanced.")
    elif unverified and len(unverified) >= n_steps:
        # Every step lacked formulas: "balanced" here would be a claim nothing supports.
        bits.append("balance ?")
        tips.append("No step could be checked — the database has no formulas for these "
                    "metabolites, so balance is unknown, not confirmed.")
    elif unverified:
        bits.append(f"balanced ✓ ({len(unverified)} unchecked)")
        tips.append(f"{len(unverified)} step(s) had no formulas and could not be checked: "
                    + ", ".join(unverified[:4]))
    else:
        bits.append("balanced ✓")
    # Name the compounds behind any "unchecked": that is the part the user can fix.
    missing = []
    for labels in (getattr(result, "unverified_reasons", None) or {}).values():
        for label in labels:
            if label not in missing:
                missing.append(label)
    if missing:
        tips.append("No formula for: " + ", ".join(missing[:4])
                    + (" …" if len(missing) > 4 else ""))
    bad = any("✗" in b for b in bits)
    unknown = any("?" in b for b in bits)
    colour = "#c5221f" if bad else ("#b06000" if unknown else "")
    return (", ".join(bits), "\n".join(tips), colour, align)


def _verdict_cell(result):
    align = Qt.AlignLeft | Qt.AlignVCenter
    try:
        from ...core import feasibility as fz
        rep = fz.assess(result,
                        diagnosis=getattr(result, "_diagnosis", None),
                        branching=getattr(result, "_branching", None),
                        mdf=getattr(result, "_mdf", None))
        return (rep.label, rep.sentence(), rep.colour, align)
    except Exception:  # noqa: BLE001 — a comparison must never fail on a verdict
        return ("—", "", style.TEXT_MUTED, align)
