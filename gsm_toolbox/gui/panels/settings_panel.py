"""The Settings tab: general app configuration.

For now it exposes the persistent **cache** — the databases and molecule-structure
images the app has downloaded — grouped by category, with the size of each file
and buttons to free space. Everything here can be re-downloaded on demand, and
keeping it lets the app work offline.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ...core import cache


# Categories shown as a single collapsed entry (total size + "Delete all") rather
# than a full per-file list — for caches with many small, uninteresting files.
_COMPACT_CATEGORIES = {"Molecule structures"}


class SettingsPanel(QWidget):
    def __init__(self):
        super().__init__()
        self._lists: dict[str, QListWidget] = {}
        self._headers: dict[str, QLabel] = {}

        title = QLabel("Settings")
        title.setStyleSheet("font-size:16px; font-weight:700; padding:4px 0;")
        self._location = QLabel()
        self._location.setStyleSheet("color:#5f6368;")
        self._location.setWordWrap(True)

        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self.refresh)
        head = QHBoxLayout()
        head.addWidget(title, 1)
        head.addWidget(refresh)

        body = QVBoxLayout()
        intro = QLabel(
            "<b>Stored cache.</b> Downloaded databases and structure images are kept "
            "on disk for offline use. Removing a file frees space; it is re-fetched "
            "if needed.")
        intro.setWordWrap(True)
        body.addLayout(head)
        body.addWidget(intro)
        body.addWidget(self._location)

        for category in cache.CATEGORIES:
            body.addWidget(self._make_category_box(category))
        body.addStretch(1)

        inner = QWidget()
        inner.setLayout(body)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(inner)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer = QVBoxLayout(self)
        outer.addWidget(scroll)

        self.refresh()

    def _make_category_box(self, category: str) -> QGroupBox:
        box = QGroupBox(category)
        v = QVBoxLayout(box)
        header = QLabel()
        header.setStyleSheet("color:#5f6368;")
        self._headers[category] = header

        if category in _COMPACT_CATEGORIES:
            # Single summary line + one "Delete all" button (no per-file list).
            buttons = QHBoxLayout()
            clear = QPushButton("Delete all")
            clear.clicked.connect(lambda _=False, c=category: self._clear_category(c))
            buttons.addWidget(header, 1)
            buttons.addStretch(1)
            buttons.addWidget(clear)
            v.addLayout(buttons)
            return box

        lst = QListWidget()
        lst.setSelectionMode(QAbstractItemView.ExtendedSelection)
        lst.setMaximumHeight(160)
        self._lists[category] = lst
        buttons = QHBoxLayout()
        del_sel = QPushButton("Delete selected")
        del_sel.clicked.connect(lambda _=False, c=category: self._delete_selected(c))
        clear = QPushButton("Clear all in this category")
        clear.clicked.connect(lambda _=False, c=category: self._clear_category(c))
        buttons.addWidget(header, 1)
        buttons.addStretch(1)
        buttons.addWidget(del_sel)
        buttons.addWidget(clear)
        v.addWidget(lst)
        v.addLayout(buttons)
        return box

    # ----- actions -----------------------------------------------------
    def refresh(self) -> None:
        self._location.setText(f"Cache location: {cache.base_dir()}")
        for category, header in self._headers.items():
            files = cache.list_category(category)
            total = sum(f["size"] for f in files)
            header.setText(f"{len(files)} file(s) · {cache.human_size(total)}")
            lst = self._lists.get(category)
            if lst is None:
                continue
            lst.clear()
            for f in files:
                item = QListWidgetItem(f"{f['name']}    ({cache.human_size(f['size'])})")
                item.setData(Qt.UserRole, f["path"])
                lst.addItem(item)

    def _delete_selected(self, category: str) -> None:
        lst = self._lists[category]
        paths = [it.data(Qt.UserRole) for it in lst.selectedItems()]
        if not paths:
            QMessageBox.information(self, "Nothing selected",
                                    "Select one or more files to delete.")
            return
        if QMessageBox.question(
                self, "Delete cached files",
                f"Delete {len(paths)} selected file(s)? They can be re-downloaded later.") \
                != QMessageBox.Yes:
            return
        freed = cache.delete_paths(paths)
        self.refresh()
        QMessageBox.information(self, "Cache cleaned",
                                f"Freed {cache.human_size(freed)}.")

    def _clear_category(self, category: str) -> None:
        files = cache.list_category(category)
        if not files:
            return
        total = sum(f["size"] for f in files)
        if QMessageBox.question(
                self, "Clear category",
                f"Delete all {len(files)} file(s) in “{category}” "
                f"({cache.human_size(total)})?") != QMessageBox.Yes:
            return
        freed = cache.clear_category(category)
        self.refresh()
        QMessageBox.information(self, "Cache cleaned",
                                f"Freed {cache.human_size(freed)}.")
