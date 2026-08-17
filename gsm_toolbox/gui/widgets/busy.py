"""A modal 'busy' dialog with an indeterminate progress bar.

Used to give feedback during blocking operations (loading a model, adding a
reaction, applying a pathway) that would otherwise freeze the interface with no
indication that anything is happening. The work runs on a background QThread so
the spinner keeps animating; the dialog closes automatically on completion.
"""

from __future__ import annotations

from typing import Any, Callable, Optional

from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtWidgets import QDialog, QLabel, QProgressBar, QVBoxLayout


class _Worker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn: Callable[[], Any]):
        super().__init__()
        self._fn = fn

    def run(self) -> None:  # noqa: D401 - QThread entry point
        try:
            result = self._fn()
        except Exception as exc:  # noqa: BLE001 - report to caller
            # Keep the traceback. Reporting only str(exc) leaves a message like
            # "No module named 'numpy._core'" with no indication of which operation
            # raised it — that one cost a great deal of guesswork to track down.
            import traceback
            detail = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
            try:
                import sys
                print(detail, file=sys.stderr, flush=True)   # → the console log
            except Exception:  # noqa: BLE001 - never fail while reporting a failure
                pass
            self.failed.emit(str(exc))
            return
        self.done.emit(result)


class BusyDialog(QDialog):
    """A small modal dialog shown while a background operation runs.

    With ``cancelable=True`` it shows a Cancel button (and re-enables the window
    close button) so the user can stop a long operation."""

    cancelled = Signal()

    def __init__(self, parent, message: str, title: str = "Please wait",
                 cancelable: bool = False):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(True)
        self._cancelable = cancelable
        self.setWindowFlag(Qt.WindowCloseButtonHint, cancelable)
        self.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        self.setMinimumWidth(340)

        layout = QVBoxLayout(self)
        self._label = QLabel(message)
        self._label.setWordWrap(True)
        layout.addWidget(self._label)

        bar = QProgressBar()
        bar.setRange(0, 0)  # indeterminate
        bar.setTextVisible(False)
        layout.addWidget(bar)

        if cancelable:
            from PySide6.QtWidgets import QDialogButtonBox
            bb = QDialogButtonBox(QDialogButtonBox.Cancel)
            bb.rejected.connect(self._request_cancel)
            layout.addWidget(bb)

    def set_message(self, message: str) -> None:
        self._label.setText(message)

    def _request_cancel(self) -> None:
        self._label.setText("Stopping…")
        self.cancelled.emit()

    def reject(self) -> None:  # noqa: D401 - Esc / window-close = cancel (if allowed)
        if self._cancelable:
            self._request_cancel()
        # otherwise ignore (dialog closes itself when the work finishes)


def run_on_gui_with_busy(parent, message: str, fn: Callable[[], Any], *,
                         title: str = "Please wait") -> Any:
    """Run GUI-thread-bound work behind a modal busy dialog.

    Some work (e.g. building a QGraphicsScene for the network map) must run on the
    GUI thread and cannot be moved to a worker. This shows the busy dialog, lets it
    paint, then runs ``fn`` synchronously and closes the dialog — so the user sees
    a "loading" popup instead of a frozen window. Returns ``fn``'s result.
    """
    from PySide6.QtWidgets import QApplication

    dlg = BusyDialog(parent, message, title=title)
    dlg.show()
    # Let the dialog and its spinner actually paint before we block the thread.
    for _ in range(3):
        QApplication.processEvents()
    try:
        return fn()
    finally:
        dlg.close()
        dlg.deleteLater()


#: returned as the result when the user cancels a cancelable run_busy.
CANCELLED = "__CANCELLED__"


def was_cancelled(result) -> bool:
    return result is CANCELLED or result == CANCELLED


class _Progress(QObject):
    """Carries progress text from the worker thread to the dialog (GUI thread)."""

    message = Signal(str)


def run_busy(parent, message: str, fn: Callable[[], Any], *,
             title: str = "Please wait",
             after: Optional[Callable[[Any], None]] = None,
             after_message: Optional[str] = None,
             cancelable: bool = False,
             progress: bool = False) -> tuple[bool, Any]:
    """Run ``fn`` on a background thread while showing a modal busy dialog.

    Returns ``(ok, result)``: ``ok`` is False and ``result`` is the error message
    string if ``fn`` raised; otherwise ``result`` is the function's return value.
    Blocks (via a modal event loop) until the work finishes, but keeps the UI
    responsive and the spinner animated.

    ``after`` (if given) is a main-thread callback run with the worker's result
    *while the dialog is still visible*, before it closes. Use it for finalisation
    that must run on the GUI thread (e.g. refreshing widgets) so the "loading"
    dialog stays up until the app is genuinely ready — the window won't be left
    frozen with the popup already gone. ``after_message`` updates the dialog text
    before that step runs.
    """
    dlg = BusyDialog(parent, message, title=title, cancelable=cancelable)
    # With progress=True, `fn` is handed a thread-safe reporter it can call to update the
    # dialog text — the worker cannot touch widgets directly, so it goes through a signal.
    reporter = None
    if progress:
        reporter = _Progress()
        reporter.message.connect(dlg.set_message)
        _emit = reporter.message.emit
        _inner = fn

        def fn():                      # noqa: F811 — deliberately rebinding the callable
            return _inner(_emit)

    worker = _Worker(fn)
    outcome: dict = {"ok": True, "result": None, "cancelled": False}

    def _cancel():
        outcome["cancelled"] = True
        outcome["ok"] = False
        outcome["result"] = CANCELLED
        if worker.isRunning():
            worker.terminate()      # the cancelable ops are read-only / build-new
        dlg.accept()
    dlg.cancelled.connect(_cancel)

    def _done(result):
        if outcome["cancelled"]:
            return
        outcome["ok"] = True
        outcome["result"] = result
        if after is not None:
            if after_message:
                dlg.set_message(after_message)
                dlg.repaint()   # show the new message before the main thread blocks
            try:
                after(result)
            except Exception as exc:  # noqa: BLE001 - surface but still close
                outcome["ok"] = False
                outcome["result"] = exc
        dlg.accept()

    def _failed(msg):
        if outcome["cancelled"]:
            return
        outcome["ok"] = False
        outcome["result"] = msg
        dlg.accept()

    worker.done.connect(_done)
    worker.failed.connect(_failed)
    worker.start()
    dlg.exec()
    worker.wait(3000)
    return outcome["ok"], outcome["result"]
