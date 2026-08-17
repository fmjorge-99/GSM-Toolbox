"""A read-only 'Check for updates' dialog.

Reports the version of the toolbox, each key dependency, and each downloadable data
asset, flags anything out of date, and — crucially — only *suggests* the manual command
to update it. Nothing here installs or mutates the environment: the user copies the
command and runs it themselves. Checking PyPI happens off the UI thread.
"""
from __future__ import annotations

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QHeaderView, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

from ...core import updater
from .. import style

_STATUS_COLOR = {
    updater.OUTDATED: "#E8710A",
    updater.MISSING: "#C5221F",
    updater.OK: "#188038",
    updater.UNKNOWN: style.TEXT_MUTED,
    updater.INFO: style.TEXT_MUTED,
}
_STATUS_LABEL = {
    updater.OUTDATED: "Update available",
    updater.MISSING: "Not installed",
    updater.OK: "Up to date",
    updater.UNKNOWN: "Unknown",
    updater.INFO: "Info",
}


class _CheckThread(QThread):
    done = Signal(object)

    def run(self):
        try:
            report = updater.build_report(online=True)
        except Exception as exc:  # noqa: BLE001
            report = exc
        self.done.emit(report)


class UpdateDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Check for updates")
        self.resize(820, 520)
        v = QVBoxLayout(self)

        self._intro = QLabel(
            "Current versions of the toolbox, its dependencies, and its data. This only "
            "<b>checks and suggests</b> — it never installs anything. Copy a suggested "
            "command into a terminal to update a component yourself.")
        self._intro.setWordWrap(True)
        v.addWidget(self._intro)

        self._status = QLabel("Checking for updates online…")
        self._status.setStyleSheet(f"color: {style.TEXT_MUTED};")
        v.addWidget(self._status)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["Component", "Installed", "Latest", "Status", "Suggested command / action"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(4, QHeaderView.Stretch)
        v.addWidget(self.table, 1)

        btns = QHBoxLayout()
        self.copy_btn = QPushButton("Copy suggested commands")
        self.copy_btn.setToolTip("Copy the update commands for every out-of-date "
                                 "component to the clipboard.")
        self.copy_btn.clicked.connect(self._copy_commands)
        self.copy_btn.setEnabled(False)
        self.recheck_btn = QPushButton("Re-check")
        self.recheck_btn.clicked.connect(self._start)
        close_btn = QPushButton("Close")
        close_btn.setObjectName("primary")
        close_btn.clicked.connect(self.accept)
        btns.addWidget(self.copy_btn)
        btns.addWidget(self.recheck_btn)
        btns.addStretch(1)
        btns.addWidget(close_btn)
        v.addLayout(btns)

        self._report = None
        self._thread = None
        self._start()

    def _start(self):
        self._status.setText("Checking for updates online…")
        self.recheck_btn.setEnabled(False)
        self.copy_btn.setEnabled(False)
        self._thread = _CheckThread(self)
        self._thread.done.connect(self._show_report)
        self._thread.start()

    def _show_report(self, report):
        self.recheck_btn.setEnabled(True)
        if isinstance(report, Exception):
            self._status.setText(f"Could not complete the check: {report}")
            return
        self._report = report
        self.table.setRowCount(0)
        for c in report.components:
            row = self.table.rowCount()
            self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(c.name))
            self.table.setItem(row, 1, QTableWidgetItem(c.current or "—"))
            self.table.setItem(row, 2, QTableWidgetItem(c.latest or "—"))
            from PySide6.QtGui import QBrush, QColor
            st = QTableWidgetItem(_STATUS_LABEL.get(c.status, c.status))
            st.setForeground(QBrush(QColor(_STATUS_COLOR.get(c.status, style.TEXT_MUTED))))
            self.table.setItem(row, 3, st)
            action = c.suggestion or (c.detail or "")
            act_item = QTableWidgetItem(action)
            if c.detail and c.suggestion:
                act_item.setToolTip(c.detail)
            self.table.setItem(row, 4, act_item)

        outdated, missing = report.n_outdated, report.n_missing
        if not report.checked_online:
            msg = ("Offline — showing installed versions only; could not check for newer "
                   "releases.")
        elif outdated or missing:
            bits = []
            if outdated:
                bits.append(f"{outdated} update(s) available")
            if missing:
                bits.append(f"{missing} required component(s) missing")
            msg = " · ".join(bits) + ". Suggested commands are listed; run them yourself."
        else:
            msg = "Everything tracked is up to date."
        self._status.setText(msg)
        self.copy_btn.setEnabled(bool(self._commands()))

    def _commands(self):
        if not self._report:
            return []
        return [c.suggestion for c in self._report.components
                if c.status in (updater.OUTDATED, updater.MISSING) and c.suggestion]

    def _copy_commands(self):
        from PySide6.QtWidgets import QApplication
        cmds = self._commands()
        if cmds:
            QApplication.clipboard().setText("\n".join(cmds))
            self._status.setText(f"Copied {len(cmds)} command(s) to the clipboard.")
