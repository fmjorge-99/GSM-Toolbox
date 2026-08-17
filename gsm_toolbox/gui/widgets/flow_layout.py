"""A wrapping ("flow") layout for toolbars of controls.

Why this exists: a long single-row ``QHBoxLayout`` of controls reports a minimum
width equal to the SUM of its children's minimums. A QTabWidget takes the maximum
over *all* its pages, so one wide control row forced the whole main window to a
minimum width larger than the screen — and because Qt always honours minimumSize
over maximumSize, the window could not stay maximized and was pushed off-screen on
every relayout.

A flow layout instead wraps its items onto as many rows as needed, so its minimum
width is only that of the widest single item. Controls therefore keep their natural
size (no cropped text) and the window is free to be any width.

Adapted from the standard Qt "Flow Layout" pattern.
"""

from __future__ import annotations

from PySide6.QtCore import QMargins, QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import QLayout, QSizePolicy, QWidget


class FlowLayout(QLayout):
    def __init__(self, parent: QWidget | None = None, margin: int = 0,
                 h_spacing: int = 6, v_spacing: int = 4):
        super().__init__(parent)
        self._items: list = []
        self._h = h_spacing
        self._v = v_spacing
        self.setContentsMargins(QMargins(margin, margin, margin, margin))

    # -- QLayout API ---------------------------------------------------
    def addItem(self, item):  # noqa: N802 - Qt override
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index):  # noqa: N802 - Qt override
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index):  # noqa: N802 - Qt override
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self):  # noqa: N802 - Qt override
        return Qt.Orientations(0)

    def hasHeightForWidth(self) -> bool:  # noqa: N802 - Qt override
        return True

    def heightForWidth(self, width: int) -> int:  # noqa: N802 - Qt override
        return self._layout(QRect(0, 0, width, 0), test_only=True)

    def setGeometry(self, rect):  # noqa: N802 - Qt override
        super().setGeometry(rect)
        self._layout(rect, test_only=False)

    def sizeHint(self) -> QSize:  # noqa: N802 - Qt override
        return self.minimumSize()

    def minimumSize(self) -> QSize:  # noqa: N802 - Qt override
        # Only as wide as the widest single item — this is what frees the window.
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        return size + QSize(m.left() + m.right(), m.top() + m.bottom())

    # -- internals -----------------------------------------------------
    def _layout(self, rect: QRect, test_only: bool) -> int:
        m = self.contentsMargins()
        eff = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        x, y, line_h = eff.x(), eff.y(), 0
        for item in self._items:
            hint = item.sizeHint()
            next_x = x + hint.width() + self._h
            if next_x - self._h > eff.right() and line_h > 0:
                x = eff.x()
                y = y + line_h + self._v
                next_x = x + hint.width() + self._h
                line_h = 0
            if not test_only:
                item.setGeometry(QRect(QPoint(x, y), hint))
            x = next_x
            line_h = max(line_h, hint.height())
        return y + line_h - rect.y() + m.bottom()
