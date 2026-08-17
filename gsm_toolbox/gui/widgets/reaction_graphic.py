"""Render a reaction graphically: substrate structures = product structures.

Main metabolites are drawn as their 2-D structure image with the name beneath;
ubiquitous currency metabolites (ATP, NADPH, H2O, ferredoxin…) are shown as a
short text label on/near the arrow rather than a structure, to keep the picture
readable — the same convention used in textbook pathway figures.
"""

from __future__ import annotations

from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ...core import network_graph as ng
from .structure_fetcher import StructureFetcher, metabolite_structure_hints


def _is_currency(met) -> bool:
    base = met.id.rsplit("_", 1)[0] if "_" in met.id else met.id
    return base.lower() in ng.CURRENCY_BASES


class _MetImage(QVBoxLayout):
    """A structure-image placeholder + name, filled asynchronously."""

    def __init__(self, met, coeff: float, size: int, fetchers: list):
        super().__init__()
        self.setAlignment(Qt.AlignCenter)
        self.img = QLabel("…")
        self.img.setAlignment(Qt.AlignCenter)
        self.img.setFixedSize(size, size)
        self.img.setStyleSheet("color:#9aa0a6;")
        label = ng.short_metabolite_name(met.id, met.name or "")
        if abs(coeff) != 1:
            label = f"{abs(coeff):g} × {label}"
        cap = QLabel(label.replace("-", "‑"))   # non-breaking hyphen: don't split "4-coumarate"
        cap.setToolTip(label)
        cap.setAlignment(Qt.AlignHCenter | Qt.AlignTop)
        cap.setWordWrap(True)
        cap.setMaximumWidth(size + 90)
        # Reserve room for a wrapped (multi-word) name so the second line is never
        # clipped — e.g. "prenyl diphosphate" must not show as just "prenyl".
        cap.setMinimumHeight(46)
        cap.setStyleSheet("font-size:11px;")
        self.addWidget(self.img, alignment=Qt.AlignHCenter)
        self.addWidget(cap, alignment=Qt.AlignHCenter)

        fetcher = StructureFetcher.for_metabolite(met, size=size)
        fetcher.fetched.connect(self._on_fetched)
        fetchers.append(fetcher)
        fetcher.start()

    def _on_fetched(self, _tag: str, data: bytes) -> None:
        if not data:
            # The name is already shown beneath; keep the tile a subtle placeholder.
            self.img.setText("—")
            self.img.setStyleSheet("color:#c0c6cc; border:1px solid #E4E8EC; border-radius:6px;")
            return
        pix = QPixmap()
        if pix.loadFromData(data):
            self.img.setPixmap(pix.scaled(self.img.width(), self.img.height(),
                                          Qt.KeepAspectRatio, Qt.SmoothTransformation))


class ReactionGraphicWidget(QScrollArea):
    """A horizontally-scrollable graphical depiction of a single reaction."""

    def __init__(self, reaction, parent=None, img_size: int = 120):
        super().__init__(parent)
        self._fetchers: List[StructureFetcher] = []
        self.setWidgetResizable(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setFrameShape(QScrollArea.NoFrame)
        self.setMinimumHeight(img_size + 66)

        inner = QWidget()
        row = QHBoxLayout(inner)
        row.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        reactants = [(m, c) for m, c in reaction.metabolites.items() if c < 0]
        products = [(m, c) for m, c in reaction.metabolites.items() if c > 0]

        self._add_side(row, reactants, img_size)
        arrow = QLabel("⇌" if reaction.reversibility else "→")
        arrow.setStyleSheet("font-size:26px; font-weight:bold; padding:0 10px;")
        row.addWidget(arrow, alignment=Qt.AlignVCenter)
        self._add_side(row, products, img_size)
        row.addStretch(1)
        self.setWidget(inner)

    def _add_side(self, row, mets, size) -> None:
        currency = [(m, c) for m, c in mets if _is_currency(m)]
        main = [(m, c) for m, c in mets if not _is_currency(m)]
        if not main and currency:      # all-currency side: show as text
            main, currency = currency, []
        first = True
        for met, coeff in main:
            if not first:
                plus = QLabel("+")
                plus.setStyleSheet("font-size:20px; padding:0 6px;")
                row.addWidget(plus, alignment=Qt.AlignVCenter)
            row.addLayout(_MetImage(met, coeff, size, self._fetchers))
            first = False
        if currency:
            names = " + ".join(
                ng.short_metabolite_name(m.id, m.name or "") for m, _ in currency)
            cur = QLabel(f"(+ {names})")
            cur.setStyleSheet("color:#5f6368; font-size:11px; padding:0 6px;")
            cur.setWordWrap(True)
            cur.setMaximumWidth(140)
            row.addWidget(cur, alignment=Qt.AlignVCenter)

    def closeEvent(self, event):  # noqa: N802 - stop threads cleanly
        for f in self._fetchers:
            f.wait(50)
        super().closeEvent(event)
