"""Status-bar job indicator + a pop-up listing all running jobs.

Shows a specific description and an estimated-progress bar for the running job;
when several jobs run at once it shows "Multiple analyses running (N)" and the
expand arrow opens a dialog detailing each one.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


def _fmt(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return f"{seconds}s"
    return f"{seconds // 60}m {seconds % 60:02d}s"


class JobsDialog(QDialog):
    """Non-modal list of all running jobs with description + estimated progress.

    Rows are kept **stable** (created once per job, then updated in place) so that
    the cancel button never gets destroyed out from under a click — rebuilding the
    whole list on every progress tick used to eat clicks and require several taps.
    """

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Running processes")
        self.resize(440, 260)
        self._manager = manager
        self._rows: dict[int, dict] = {}

        self._layout = QVBoxLayout(self)
        self._empty = QLabel("No analyses are currently running.")
        self._empty.setStyleSheet("color:#5f6368; padding:8px;")
        self._layout.addWidget(self._empty)
        self._rows_box = QVBoxLayout()
        self._layout.addLayout(self._rows_box)
        self._layout.addStretch(1)
        self.refresh()

    def _make_row(self, jid: int) -> dict:
        box = QWidget()
        v = QVBoxLayout(box)
        v.setContentsMargins(6, 6, 6, 6)
        head = QHBoxLayout()
        title = QLabel()
        title.setWordWrap(True)
        head.addWidget(title, 1)
        cancel = QToolButton()
        cancel.setText("✕")
        cancel.setAutoRaise(True)
        cancel.setStyleSheet("QToolButton{color:#c0392b; font-weight:bold;}")
        # Connected once and never recreated, so a single click always lands.
        cancel.clicked.connect(lambda _=False, i=jid: self._cancel(i))
        head.addWidget(cancel)
        v.addLayout(head)
        bar = QProgressBar()
        bar.setRange(0, 100)
        v.addWidget(bar)
        sub = QLabel()
        sub.setStyleSheet("color:#5f6368;")
        v.addWidget(sub)
        divider = QFrame()
        divider.setFrameShape(QFrame.HLine)
        divider.setStyleSheet("color:#EEF0F3;")
        self._rows_box.addWidget(box)
        self._rows_box.addWidget(divider)
        return {"box": box, "divider": divider, "title": title, "cancel": cancel,
                "bar": bar, "sub": sub}

    def _cancel(self, jid: int) -> None:
        self._manager.cancel(jid)
        self.refresh()

    def refresh(self) -> None:
        jobs = {job["id"]: job for job in self._manager.snapshot()}
        self._empty.setVisible(not jobs)

        # Remove rows for jobs that are gone.
        for jid in [j for j in self._rows if j not in jobs]:
            row = self._rows.pop(jid)
            row["box"].deleteLater()
            row["divider"].deleteLater()

        for jid, job in jobs.items():
            row = self._rows.get(jid)
            if row is None:
                row = self._rows[jid] = self._make_row(jid)
            queued = job["state"] == "queued"
            cancelling = job["state"] == "cancelling"
            row["title"].setText(f"<b>{job['title']}</b>")
            row["cancel"].setToolTip("Remove from queue" if queued else "Cancel this analysis")
            row["cancel"].setEnabled(not cancelling)
            row["bar"].setValue(0 if queued else int(job["fraction"] * 100))
            if queued:
                row["sub"].setText("queued — will start when a slot is free")
            elif cancelling:
                row["sub"].setText("cancelling…")
            else:
                row["sub"].setText(f"running · elapsed {_fmt(job['elapsed'])} · "
                                   f"about {_fmt(job['eta'])} left (estimated)")


class JobStatusBar(QWidget):
    """A status-bar widget: description (left) + progress + expand arrow (right)."""

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self._manager = manager
        self._idle_text = "Ready. Open a model to begin (File ▸ Open Model)."
        self._dialog: JobsDialog | None = None

        self._label = QLabel(self._idle_text)
        self._label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        # A status message must NEVER be able to resize the window. A plain QLabel
        # reports its full text width as its size hint, so a long message ("Saved
        # strategy 'Round 2 - …'", "FBA complete: objective = …") pushed the status
        # bar's minimum past the screen width — and because Qt honours minimumSize
        # over the maximized state, the window silently dropped out of Maximized.
        # That is why it happened after almost any analysis: they all post a message.
        # Ignored policy = "my text does not constrain my width"; the text is elided
        # to whatever room the bar actually has (see _apply_elided_text).
        self._label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self._label.setMinimumWidth(0)
        self._full_text = self._idle_text
        self._bar = QProgressBar()
        self._bar.setMaximumWidth(180)
        self._bar.setRange(0, 100)
        self._bar.hide()
        self._expand = QToolButton()
        self._expand.setText("▸")
        self._expand.setToolTip("Show all running processes")
        self._expand.setAutoRaise(True)
        self._expand.hide()
        self._expand.clicked.connect(self._toggle_dialog)

        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self._label, 1)
        lay.addWidget(self._bar)
        lay.addWidget(self._expand)

        self._manager.updated.connect(self.refresh)

    # The main window calls setText(...) to post idle/completion messages,
    # so existing call sites keep working unchanged.
    def setText(self, text: str) -> None:  # noqa: N802
        self._idle_text = text
        self.refresh()

    def message(self, text: str) -> None:
        self.setText(text)

    def _set_message(self, text: str) -> None:
        """Show `text`, elided to the space available, with the full text on hover."""
        self._full_text = text or ""
        self._label.setToolTip(self._full_text)
        self._apply_elided_text()

    def _apply_elided_text(self) -> None:
        avail = max(0, self._label.width())
        if avail <= 0:
            self._label.setText(self._full_text)
            return
        fm = self._label.fontMetrics()
        self._label.setText(fm.elidedText(self._full_text, Qt.ElideRight, avail))

    def resizeEvent(self, event):  # noqa: N802 - Qt override
        # Re-elide to the new width: the label is Ignored-policy, so it is given
        # whatever room the bar has rather than demanding room for its text.
        super().resizeEvent(event)
        self._apply_elided_text()

    def refresh(self) -> None:
        jobs = self._manager.snapshot()
        n = len(jobs)
        running = [j for j in jobs if j["state"] in ("running", "cancelling")]
        queued = [j for j in jobs if j["state"] == "queued"]
        if n == 0:
            self._set_message(self._idle_text)
            self._bar.hide()
            self._expand.hide()
        elif n == 1 and running:
            job = running[0]
            self._set_message(job["title"])
            self._bar.show()
            self._bar.setValue(int(job["fraction"] * 100))
            self._bar.setToolTip(f"about {_fmt(job['eta'])} left (estimated)")
            self._expand.show()
        else:
            q = f", {len(queued)} queued" if queued else ""
            self._set_message(f"Multiple analyses running ({len(running)}{q})")
            self._bar.show()
            self._bar.setValue(int(min((j["fraction"] for j in running), default=0.0) * 100)
                               if running else 0)
            self._bar.setToolTip("Click the arrow for per-process details and to cancel")
            self._expand.show()
        if self._dialog is not None and self._dialog.isVisible():
            self._dialog.refresh()

    def _toggle_dialog(self) -> None:
        if self._dialog is None:
            self._dialog = JobsDialog(self._manager, self.window())
        if self._dialog.isVisible():
            self._dialog.hide()
            self._expand.setText("▸")
        else:
            self._dialog.refresh()
            self._dialog.show()
            self._dialog.raise_()
            self._expand.setText("▾")
