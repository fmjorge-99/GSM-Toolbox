"""A push button whose label word-wraps to the available width and grows tall
enough to show every line — so text is never cropped, and the number of lines
adapts as the panel is resized (#T3)."""

from __future__ import annotations

from PySide6.QtWidgets import QPushButton, QSizePolicy


class WrapButton(QPushButton):
    def __init__(self, text: str = "", parent=None):
        super().__init__(parent)
        self._full = text
        super().setText(text)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        self._last_w = -1

    def setText(self, text: str) -> None:  # noqa: N802 - Qt override
        self._full = text
        super().setText(text)
        self._last_w = -1
        self._rewrap()

    def resizeEvent(self, event):  # noqa: N802 - Qt override
        super().resizeEvent(event)
        self._rewrap()

    def _rewrap(self) -> None:
        avail = self.width() - 22        # room for horizontal padding
        if avail < 30 or avail == self._last_w:
            return
        self._last_w = avail
        fm = self.fontMetrics()
        words = self._full.split()
        lines, cur = [], (words[0] if words else self._full)
        for w in words[1:]:
            candidate = f"{cur} {w}"
            if fm.horizontalAdvance(candidate) <= avail:
                cur = candidate
            else:
                lines.append(cur)
                cur = w
        lines.append(cur)
        super().setText("\n".join(lines))
        needed = fm.height() * len(lines) + 18
        if self.minimumHeight() != needed:
            self.setMinimumHeight(needed)
