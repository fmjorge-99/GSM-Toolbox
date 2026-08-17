"""Helpers to keep pop-up windows usable — never larger than the screen —
plus shared search-combo and file-dialog helpers that fix the "typed text is
appended instead of replacing/filtering" class of bugs (#2 / #7)."""

from __future__ import annotations

import os

from PySide6.QtCore import QEvent, QObject, Qt
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QCompleter,
    QFileDialog,
    QWidget,
)


class _SelectAllOnFocus(QObject):
    """Event filter that selects a line edit's whole text on focus-in, so the
    user's first keystroke REPLACES the current value instead of appending to it."""

    def eventFilter(self, obj, event):  # noqa: N802 - Qt override
        if event.type() == QEvent.FocusIn:
            from PySide6.QtCore import QTimer
            QTimer.singleShot(0, obj.selectAll)   # after Qt's own focus handling
        return False


def configure_search_combo(combo: QComboBox) -> None:
    """Turn an editable combo into a forgiving search box: a contains-match popup
    completer over its items, and select-all-on-focus so typing replaces the text
    (fixes the 'EX_ac_eEX_etoh' concatenation and prefix-jump, #2)."""
    combo.setEditable(True)
    combo.setInsertPolicy(QComboBox.NoInsert)
    combo.setMaxVisibleItems(20)
    completer = QCompleter(combo.model(), combo)
    completer.setCaseSensitivity(Qt.CaseInsensitive)
    completer.setFilterMode(Qt.MatchContains)
    completer.setCompletionMode(QCompleter.PopupCompletion)
    combo.setCompleter(completer)
    line = combo.lineEdit()
    if line is not None:
        # Keep a reference on the combo so the filter isn't garbage-collected.
        combo._focus_filter = _SelectAllOnFocus(combo)
        line.installEventFilter(combo._focus_filter)


def _default_documents_dir() -> str:
    for name in ("Documents", "Desktop"):
        d = os.path.join(os.path.expanduser("~"), name)
        if os.path.isdir(d):
            return d
    return os.path.expanduser("~")


def choose_save_path(parent, caption: str, default_name: str, name_filter: str,
                     start_dir: str = "") -> str:
    """Save-file dialog that separates the starting directory from the default file
    name (never passing the name as the ``dir`` argument), so a typed absolute path
    can't be concatenated onto the default name (#7). Returns "" if cancelled."""
    dlg = QFileDialog(parent, caption)
    dlg.setAcceptMode(QFileDialog.AcceptSave)
    dlg.setFileMode(QFileDialog.AnyFile)
    if name_filter:
        dlg.setNameFilter(name_filter)
    dlg.setDirectory(start_dir or _default_documents_dir())
    if default_name:
        dlg.selectFile(default_name)
    if dlg.exec() == QFileDialog.Accepted:
        files = dlg.selectedFiles()
        return files[0] if files else ""
    return ""


def choose_open_path(parent, caption: str, name_filter: str, start_dir: str = "") -> str:
    """Open-file dialog counterpart to :func:`choose_save_path`."""
    dlg = QFileDialog(parent, caption)
    dlg.setAcceptMode(QFileDialog.AcceptOpen)
    dlg.setFileMode(QFileDialog.ExistingFile)
    if name_filter:
        dlg.setNameFilter(name_filter)
    dlg.setDirectory(start_dir or _default_documents_dir())
    if dlg.exec() == QFileDialog.Accepted:
        files = dlg.selectedFiles()
        return files[0] if files else ""
    return ""


def clamp_to_screen(widget: QWidget, frac: float = 0.9) -> None:
    """Cap a window's size to a fraction of the available screen and re-centre it,
    so a content-sized dialog can never open bigger than the display (which would
    make it impossible to move or close)."""
    screen = widget.screen() or QApplication.primaryScreen()
    if screen is None:
        return
    geo = screen.availableGeometry()
    max_w, max_h = int(geo.width() * frac), int(geo.height() * frac)
    widget.setMaximumSize(max_w, max_h)
    w = min(widget.width(), max_w)
    h = min(widget.height(), max_h)
    widget.resize(w, h)
    widget.move(geo.x() + (geo.width() - w) // 2, geo.y() + (geo.height() - h) // 2)
