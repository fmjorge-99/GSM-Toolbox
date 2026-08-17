"""The Summary dashboard tab: model counts, objective, last growth rate, edits.

Wrapped in a scroll area so the content stays fully readable even when the
Information panel is docked as a short bottom strip (where a plain form layout
would otherwise clip its rows).
"""

from __future__ import annotations

from typing import Optional

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QGroupBox,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ...core import io_models


class SummaryView(QWidget):
    def __init__(self):
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QScrollArea.NoFrame)
        outer.addWidget(self._scroll)

        self._content = QWidget()
        self._layout = QVBoxLayout(self._content)
        self._layout.setAlignment(Qt.AlignTop)
        self._scroll.setWidget(self._content)

        self._show_placeholder()

    # -- internals -----------------------------------------------------
    def _clear(self) -> None:
        while self._layout.count():
            item = self._layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

    def _show_placeholder(self) -> None:
        ph = QLabel("Open a model to see its summary.")
        ph.setStyleSheet("color: #5f6368; padding: 12px;")
        ph.setWordWrap(True)
        self._layout.addWidget(ph)

    @staticmethod
    def _row(label: str, value) -> QLabel:
        """One label:value row as a single word-wrapped rich-text line.

        Rendered as one QLabel (not a two-column form) so the value can never be
        clipped to zero width when the Information panel is docked as a narrow or
        short strip (fixes blank values, #B1)."""
        v = QLabel(f"<b>{label}:</b> {value}")
        v.setTextFormat(Qt.RichText)
        v.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.setWordWrap(True)
        return v

    def _box(self, title: str, rows) -> QGroupBox:
        box = QGroupBox(title)
        v = QVBoxLayout(box)
        v.setSpacing(2)
        for label, value in rows:
            v.addWidget(self._row(label, value))
        return box

    # -- public API ----------------------------------------------------
    def update_summary(self, model: Optional[cobra.Model], diff: dict | None = None,
                       last_growth: Optional[float] = None) -> None:
        self._clear()

        if model is None:
            self._show_placeholder()
            return

        summary = io_models.summarize(model)
        self._layout.addWidget(self._box("Model", summary.as_rows()))

        if last_growth is not None:
            self._layout.addWidget(
                self._box("Last simulation", [("Objective value", f"{last_growth:.6g}")]))

        if diff:
            added = diff.get("added_reactions", [])
            removed = diff.get("removed_reactions", [])
            changed = diff.get("changed_bounds", [])
            if added or removed or changed:
                self._layout.addWidget(self._box("Changes vs original model", [
                    ("Reactions added", len(added)),
                    ("Reactions removed", len(removed)),
                    ("Bounds changed", len(changed)),
                ]))
