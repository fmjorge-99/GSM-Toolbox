"""A section that can be folded away, for dialogs whose explanations outgrow the screen.

Settings panels tend to accumulate: each option earns a paragraph explaining what it costs
and when to use it, and after three or four options the dialog is taller than the display
— at which point the OK button is off-screen and the window cannot be dismissed.

Folding the prose away solves that without deleting it. The header always shows the
control itself (a checkbox, a spin box) plus a one-line summary, so the dialog is usable
fully collapsed; the arrow reveals the detail for anyone who wants it.
"""
from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QSizePolicy, QToolButton, QVBoxLayout, QWidget)

from .. import style


class CollapsibleSection(QFrame):
    """A titled section whose body can be shown or hidden.

    Put the always-relevant control in :attr:`header_row` and the explanation in
    :attr:`body`; the body starts collapsed unless ``expanded=True``.
    """

    def __init__(self, title: str, parent: Optional[QWidget] = None, *,
                 expanded: bool = False, summary: str = ""):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        self.setObjectName("collapsibleSection")

        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 6, 8, 6)
        outer.setSpacing(4)

        self._toggle = QToolButton()
        self._toggle.setStyleSheet("QToolButton { border: none; font-weight: 600; }")
        self._toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._toggle.setArrowType(Qt.DownArrow if expanded else Qt.RightArrow)
        self._toggle.setText(title)
        self._toggle.setCheckable(True)
        self._toggle.setChecked(expanded)
        self._toggle.setCursor(Qt.PointingHandCursor)
        self._toggle.setToolTip("Show or hide the details for this section")
        self._toggle.toggled.connect(self._on_toggled)

        title_row = QHBoxLayout()
        title_row.setContentsMargins(0, 0, 0, 0)
        title_row.addWidget(self._toggle)
        self._summary = QLabel(summary)
        self._summary.setStyleSheet(f"color:{style.TEXT_MUTED};")
        self._summary.setVisible(bool(summary) and not expanded)
        title_row.addWidget(self._summary, 1)
        outer.addLayout(title_row)

        #: Always visible — put the actual control here, not in the body.
        self.header_row = QVBoxLayout()
        self.header_row.setContentsMargins(16, 0, 0, 0)
        self.header_row.setSpacing(4)
        outer.addLayout(self.header_row)

        #: Collapsible — put the explanation and any rarely-touched controls here.
        self.body = QWidget()
        self.body.setVisible(expanded)
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(16, 4, 0, 0)
        self.body_layout.setSpacing(6)
        self.body.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        outer.addWidget(self.body)

    def _on_toggled(self, on: bool) -> None:
        self._toggle.setArrowType(Qt.DownArrow if on else Qt.RightArrow)
        self.body.setVisible(on)
        self._summary.setVisible(bool(self._summary.text()) and not on)
        # Resize the dialog to match the new content — but ONLY when nothing is
        # scrolling us. Inside a scroll area the window's size hint is unrelated to the
        # content height, so adjustSize() would shrink the dialog every time a section
        # is expanded, which is precisely backwards. There the scroll bar is the
        # correct answer to overflow and the window should simply stay put.
        if self._inside_scroll_area():
            return
        window = self.window()
        if window is not None:
            window.adjustSize()

    def _inside_scroll_area(self) -> bool:
        from PySide6.QtWidgets import QAbstractScrollArea
        node = self.parentWidget()
        while node is not None:
            if isinstance(node, QAbstractScrollArea):
                return True
            node = node.parentWidget()
        return False

    def add_widget(self, widget: QWidget, *, always_visible: bool = False) -> None:
        (self.header_row if always_visible else self.body_layout).addWidget(widget)

    def add_layout(self, layout, *, always_visible: bool = False) -> None:
        (self.header_row if always_visible else self.body_layout).addLayout(layout)

    def set_expanded(self, on: bool) -> None:
        self._toggle.setChecked(bool(on))

    def is_expanded(self) -> bool:
        return self._toggle.isChecked()
