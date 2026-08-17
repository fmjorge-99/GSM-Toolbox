"""A classic, saveable metabolic-pathway scheme.

Draws the heterologous pathway top-to-bottom: the nearest native precursor at the
top, each main carbon skeleton as its 2-D structure image with the name beneath,
joined by downward arrows labelled with the enzyme/reaction. Currency co-substrates
and co-products (ATP, NADPH, ferredoxin, CoA…) are written as short names beside
each arrow; a reaction that fuses in a second carbon backbone shows that partner as
a small structure feeding into the arrow. The image can be saved to a file.

The reaction band height is sized to the (word-wrapped) enzyme text so names are
never cropped, and the enzyme text sits clear of the structure images. Long or
multi-step pathways wrap across columns to stay within a 2:1 aspect ratio.
"""

from __future__ import annotations

import math
from typing import List

from PySide6.QtCore import QRect, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QImage, QPainter, QPen, QPixmap
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
)

_S = 120          # structure image cell (px)
_NAME_H = 62      # room for the metabolite name under each image (allow ~3 lines)
_COL_W = 470      # width of one column (image in the middle, labels either side)
_MARGIN = 28
_PAD = 10


class PathwaySchemeDialog(QDialog):
    def __init__(self, parent, nodes: List[dict], arrows: List[dict], target: str = ""):
        super().__init__(parent)
        self.setWindowTitle(f"Pathway scheme — {target}" if target else "Pathway scheme")
        self.resize(760, 820)
        self._image = self._render(nodes, arrows)

        layout = QVBoxLayout(self)
        caption = QLabel("A metabolic-pathway view of the heterologous route. Save it as an "
                         "image for figures or notebooks.")
        caption.setWordWrap(True)
        layout.addWidget(caption)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        canvas = QLabel()
        canvas.setAlignment(Qt.AlignCenter)
        canvas.setPixmap(QPixmap.fromImage(self._image))
        scroll.setWidget(canvas)
        layout.addWidget(scroll, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Close, self)
        save = QPushButton("Save image…")
        save.setObjectName("primary")
        save.clicked.connect(self._save)
        buttons.addButton(save, QDialogButtonBox.ActionRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

    # ----- rendering ---------------------------------------------------
    def _zone_w(self) -> int:
        return _COL_W // 2 - _S // 2 - 22

    def _arrow_band_height(self, arrows, fm_enz, fm_cur) -> int:
        """Uniform reaction-band height sized to the tallest label so nothing is
        cropped and images never overlap the text."""
        zone = self._zone_w()
        need = 96
        for a in arrows:
            right = fm_enz.boundingRect(QRect(0, 0, zone, 9000),
                                        Qt.TextWordWrap, a.get("enzyme", "")).height()
            if a.get("produced"):
                right += fm_cur.boundingRect(QRect(0, 0, zone, 9000),
                                             Qt.TextWordWrap, "→ " + a["produced"]).height() + 4
            left = 0
            for m in a.get("merges", []):
                left += 70 + fm_cur.boundingRect(QRect(0, 0, zone, 9000),
                                                 Qt.TextWordWrap, m.get("name", "")).height() + 6
            if a.get("consumed"):
                left += fm_cur.boundingRect(QRect(0, 0, zone, 9000),
                                            Qt.TextWordWrap, "+ " + a["consumed"]).height() + 4
            need = max(need, right + 20, left + 20)
        return min(need, 360)

    def _render(self, nodes: List[dict], arrows: List[dict]) -> QImage:
        n = max(1, len(nodes))
        enz_font = QFont(); enz_font.setPointSize(10); enz_font.setBold(True)
        cur_font = QFont(); cur_font.setPointSize(9)
        name_font = QFont(); name_font.setPointSize(10); name_font.setBold(True)
        fm_enz, fm_cur = QFontMetrics(enz_font), QFontMetrics(cur_font)

        arrow_h = self._arrow_band_height(arrows, fm_enz, fm_cur)
        row_h = _S + _NAME_H + arrow_h

        cols = 1
        while cols < n:
            rows = math.ceil(n / cols)
            w = _MARGIN * 2 + cols * _COL_W
            h = _MARGIN * 2 + rows * row_h
            if h <= 2 * w:
                break
            cols += 1
        rows = math.ceil(n / cols)
        width = _MARGIN * 2 + cols * _COL_W
        height = _MARGIN * 2 + rows * row_h

        img = QImage(width, height, QImage.Format_ARGB32)
        img.fill(QColor("#FFFFFF"))
        p = QPainter(img)
        p.setRenderHint(QPainter.Antialiasing, True)
        p.setRenderHint(QPainter.SmoothPixmapTransform, True)

        def col_row(i):
            return i // rows, i % rows

        def colx(c):
            return _MARGIN + c * _COL_W

        def node_cx(c):
            return colx(c) + _COL_W // 2

        def node_top(r):
            return _MARGIN + r * row_h

        for i, node in enumerate(nodes):
            c, r = col_row(i)
            cx, y = node_cx(c), node_top(r)
            self._draw_structure(p, node.get("img"), cx - _S // 2, y, _S, cur_font,
                                 node.get("name", ""))
            p.setPen(QColor("#202124")); p.setFont(name_font)
            p.drawText(QRect(colx(c) + 6, y + _S, _COL_W - 12, _NAME_H),
                       Qt.AlignHCenter | Qt.AlignTop | Qt.TextWordWrap, node.get("name", ""))

        for i, a in enumerate(arrows):
            if i + 1 >= n:
                break
            c0, r0 = col_row(i)
            c1, r1 = col_row(i + 1)
            if c0 == c1:
                self._varrow(p, node_cx(c0), colx(c0), node_top(r0) + _S + _NAME_H,
                             node_top(r1), a, enz_font, cur_font)
            else:
                self._elbow(p, node_cx(c0), node_top(r0) + _S + _NAME_H,
                            node_cx(c1), node_top(r1), height, a, enz_font)
        p.end()
        return img

    def _draw_structure(self, p, data, x, y, size, cur_font, name: str = "") -> None:
        pix = QPixmap()
        if data and pix.loadFromData(data):
            scaled = pix.scaled(size, size, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            p.drawPixmap(x + (size - scaled.width()) // 2,
                         y + (size - scaled.height()) // 2, scaled)
            return
        # No structure available: draw a clean labelled tile with the compound name
        # (and, for a carrier conjugate, the core + a '–CoA'-style chip) so the map
        # stays readable instead of showing an empty box.
        p.setPen(QPen(QColor("#CBD2D9"))); p.setBrush(QColor("#F7F9FB"))
        p.drawRoundedRect(QRect(x, y, size, size), 8, 8)
        p.setBrush(Qt.NoBrush)
        from ..widgets.structure_fetcher import split_cofactor
        core, tag = split_cofactor(name)
        label = (name or "?").replace("-", "‑")
        chip = ""
        if tag:
            label = (core or "?").replace("-", "‑")
            chip = f"–{tag}"
        f = QFont(cur_font); f.setPointSize(11)
        p.setFont(f); p.setPen(QColor("#5f6368"))
        p.drawText(QRect(x + 6, y + 6, size - 12, size - 28),
                   Qt.AlignCenter | Qt.TextWordWrap, label)
        if chip:
            f2 = QFont(cur_font); f2.setBold(True); p.setFont(f2)
            p.setPen(QColor("#1a73e8"))
            p.drawText(QRect(x + 6, y + size - 24, size - 12, 20),
                       Qt.AlignRight | Qt.AlignVCenter, chip)

    def _varrow(self, p, cx, colx, top, bot, a, enz_font, cur_font) -> None:
        zone = self._zone_w()
        band = bot - top
        # arrow shaft
        p.setPen(QPen(QColor("#202124"), 2))
        p.drawLine(cx, top + 6, cx, bot - 12)
        p.drawLine(cx, bot - 12, cx - 6, bot - 22)
        p.drawLine(cx, bot - 12, cx + 6, bot - 22)

        # RIGHT: enzyme name (blue) then co-products (purple)
        rx = cx + _S // 2 + 12
        fm_enz = QFontMetrics(enz_font)
        enz_h = fm_enz.boundingRect(QRect(0, 0, zone, 9000), Qt.TextWordWrap,
                                    a.get("enzyme", "")).height()
        p.setFont(enz_font); p.setPen(QColor("#1a73e8"))
        p.drawText(QRect(rx, top + 6, zone, band - 12),
                   Qt.AlignLeft | Qt.AlignTop | Qt.TextWordWrap, a.get("enzyme", ""))
        if a.get("produced"):
            p.setFont(cur_font); p.setPen(QColor("#9334e6"))
            p.drawText(QRect(rx, top + 10 + enz_h, zone, band - enz_h - 14),
                       Qt.AlignLeft | Qt.AlignTop | Qt.TextWordWrap, "→ " + a["produced"])

        # LEFT: merge-in structures (fused backbones) then co-substrates (green)
        lx0 = colx + _PAD
        ly = top + 6
        for m in a.get("merges", []):
            self._draw_structure(p, m.get("img"), cx - _S // 2 - 12 - 66, ly, 66, cur_font)
            p.setFont(cur_font); p.setPen(QColor("#188038"))
            nb = QRect(lx0, ly + 66, zone, 30)
            p.drawText(nb, Qt.AlignLeft | Qt.AlignTop | Qt.TextWordWrap, "+ " + m.get("name", ""))
            ly += 70 + QFontMetrics(cur_font).boundingRect(
                nb, Qt.TextWordWrap, "+ " + m.get("name", "")).height() + 6
        if a.get("consumed"):
            p.setFont(cur_font); p.setPen(QColor("#188038"))
            p.drawText(QRect(lx0, ly, zone, bot - ly - 4),
                       Qt.AlignLeft | Qt.AlignTop | Qt.TextWordWrap, "+ " + a["consumed"])

    def _elbow(self, p, x0, y0, x1, y1, height, a, enz_font) -> None:
        yb = height - _MARGIN
        xg = x1 - _COL_W // 2 + 14
        ymid = y1 + _S // 2
        p.setPen(QPen(QColor("#9aa0a6"), 2, Qt.DashLine))
        p.drawLine(x0, y0 + 4, x0, yb)
        p.drawLine(x0, yb, xg, yb)
        p.drawLine(xg, yb, xg, ymid)
        p.setPen(QPen(QColor("#5f6368"), 2))
        xin = x1 - _S // 2 - 2
        p.drawLine(xg, ymid, xin, ymid)
        p.drawLine(xin, ymid, xin - 10, ymid - 6)
        p.drawLine(xin, ymid, xin - 10, ymid + 6)
        p.setFont(enz_font); p.setPen(QColor("#1a73e8"))
        p.drawText(QRect(x0 + 8, yb - 40, max(80, xg - x0 - 12), 38),
                   Qt.AlignHCenter | Qt.AlignBottom | Qt.TextWordWrap, a.get("enzyme", ""))

    def _save(self) -> None:
        from ..widgets.dialog_util import choose_save_path
        path = choose_save_path(self, "Save pathway scheme", "pathway_scheme.png",
                                "PNG image (*.png);;JPEG image (*.jpg)")
        if path:
            if not path.lower().endswith((".png", ".jpg", ".jpeg")):
                path += ".png"
            self._image.save(path)
