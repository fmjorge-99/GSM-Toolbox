"""Left-dock Categories panel: define and manage colored groups of reactions.

A category groups reactions of a pathway/module the user cares about. From here
they can add the reactions currently selected in the explorer, view a category in
isolation on the map, or run an analysis restricted to it.
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ...core.categories import CategoryManager


def _color_icon(hex_color: str) -> QIcon:
    pix = QPixmap(14, 14)
    pix.fill(QColor(hex_color))
    return QIcon(pix)


class CategoriesPanel(QWidget):
    new_requested = Signal()
    delete_requested = Signal(str)
    add_selected_requested = Signal(str)
    remove_selected_requested = Signal(str)
    isolate_requested = Signal(str)
    analyze_requested = Signal(str)

    def __init__(self):
        super().__init__()
        intro = QLabel("Group reactions into categories to explore pathways and analyze subsets.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #6B7280;")

        self.list = QListWidget()
        self.list.setIconSize(QSize(14, 14))

        new_btn = QPushButton("New…")
        new_btn.setObjectName("primary")
        new_btn.clicked.connect(self.new_requested)
        del_btn = QPushButton("Delete")
        del_btn.clicked.connect(self._emit_delete)

        add_btn = QPushButton("Add selected")
        add_btn.setToolTip("Add the reactions currently selected in the explorer to this category.")
        add_btn.clicked.connect(self._emit_add)
        rem_btn = QPushButton("Remove selected")
        rem_btn.clicked.connect(self._emit_remove)

        isolate_btn = QPushButton("Isolate in map")
        isolate_btn.setToolTip("Show only this category's reactions on the network map.")
        isolate_btn.clicked.connect(self._emit_isolate)
        analyze_btn = QPushButton("Analyze subset")
        analyze_btn.setToolTip("Open the Analysis tab with this category as the scope.")
        analyze_btn.clicked.connect(self._emit_analyze)

        grid = QGridLayout()
        grid.addWidget(new_btn, 0, 0)
        grid.addWidget(del_btn, 0, 1)
        grid.addWidget(add_btn, 1, 0)
        grid.addWidget(rem_btn, 1, 1)
        grid.addWidget(isolate_btn, 2, 0)
        grid.addWidget(analyze_btn, 2, 1)

        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addWidget(self.list, 1)
        layout.addLayout(grid)

    def refresh(self, manager: CategoryManager) -> None:
        current = self.selected_category_name()
        self.list.clear()
        for cat in manager.all():
            item = QListWidgetItem(_color_icon(cat.color), f"{cat.name}  ({len(cat.reaction_ids)})")
            item.setData(Qt.UserRole, cat.name)
            self.list.addItem(item)
        if current:
            self.select_category(current)

    def select_category(self, name: str) -> None:
        for i in range(self.list.count()):
            if self.list.item(i).data(Qt.UserRole) == name:
                self.list.setCurrentRow(i)
                return

    def selected_category_name(self) -> Optional[str]:
        item = self.list.currentItem()
        return None if item is None else item.data(Qt.UserRole)

    # emit helpers (no-op if nothing selected) --------------------------
    def _with_selection(self, signal):
        name = self.selected_category_name()
        if name:
            signal.emit(name)

    def _emit_delete(self):
        self._with_selection(self.delete_requested)

    def _emit_add(self):
        self._with_selection(self.add_selected_requested)

    def _emit_remove(self):
        self._with_selection(self.remove_selected_requested)

    def _emit_isolate(self):
        self._with_selection(self.isolate_requested)

    def _emit_analyze(self):
        self._with_selection(self.analyze_requested)
