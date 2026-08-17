"""One stored dynamic run, opened in its own window.

Results used to share the Dynamic Analysis panel with the settings that produced them,
which left both cramped — a table collapsed to a row or two above parameter rows squeezed
to illegibility. Giving results a window of their own lets the settings use the whole
panel and lets a table be a table.

The window is non-modal so several runs can be open side by side, which is the comparison
the run tabs exist to support.
"""
from __future__ import annotations

import re
from typing import Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QLabel, QPushButton, QTextEdit, QVBoxLayout)

from .. import style
from ..views.results_view import ResultsView


class RunResultsDialog(QDialog):
    """The table, the commentary, and a way straight to the plot."""

    #: Ask the owner to plot this run, by name.
    plot_requested = Signal(str)

    def __init__(self, parent, name: str, record: dict):
        super().__init__(parent)
        self._name = name
        self.setWindowTitle(f"{name} — {record.get('title', 'Results')}")
        self.resize(1000, 620)
        # Non-modal: two runs open at once is the point of storing them separately.
        self.setModal(False)
        self.setAttribute(Qt.WA_DeleteOnClose)

        outer = QVBoxLayout(self)

        warning = record.get("warning") or ""
        if warning:
            note = QLabel(warning)
            note.setTextFormat(Qt.RichText)
            note.setWordWrap(True)
            note.setToolTip(re.sub(r"<[^>]+>", "", warning))
            outer.addWidget(note)

        frame = record.get("frame")
        self.results = ResultsView()
        self.results.show_dataframe(frame, record.get("title", name))
        outer.addWidget(self.results, 1)

        commentary = record.get("commentary") or ""
        if commentary:
            text = QTextEdit()
            text.setReadOnly(True)
            text.setPlainText(commentary)
            text.setMaximumHeight(160)
            outer.addWidget(text, 0)

        rows = 0 if frame is None else len(frame)
        summary = QLabel(f"{rows} row(s)")
        summary.setStyleSheet(f"color:{style.TEXT_MUTED};")

        buttons = QHBoxLayout()
        plot = QPushButton("Plot Results")
        plot.setObjectName("primary")
        plot.clicked.connect(lambda: self.plot_requested.emit(self._name))
        buttons.addWidget(plot)
        buttons.addWidget(summary)
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        buttons.addWidget(close)
        outer.addLayout(buttons)

        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
