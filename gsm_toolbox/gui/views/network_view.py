"""Interactive metabolic network map (QGraphicsView).

A navigation tool: click a node to inspect it (and sync the lists), right-click a
node for actions, focus on a reaction's neighborhood or a whole category, drag
nodes around, overlay FBA fluxes, and export the map as a PNG.

Node iconography encodes biology at a glance:
  * metabolites  — ball-and-stick molecules,
  * reversible reactions   — hexagon with a double-headed (↔) arrow,
  * irreversible reactions — hexagon with a single (→) arrow,
  * transport reactions    — a membrane badge marked "T",
  * exchange reactions     — a boundary badge marked "E".

A side "Details" panel shows the selected item's properties: a metabolite's
formula, charge, identifiers and 2-D structure, or a reaction's equation, EC
number, direction and current flux (after an analysis).
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import cobra
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QImage,
    QPainter,
    QPen,
    QPixmap,
    QPolygonF,
)
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGraphicsItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from ...core import network_graph
from ...core.databases import reaction_ec_numbers
from .. import style
from ..widgets.structure_fetcher import StructureFetcher

_SCALE = 620.0
_DATA_LABEL = 0
_DATA_KIND = 1
_LABEL_FONT = QFont("Segoe UI", 8)
_LABEL_FONT.setStyleStrategy(QFont.PreferAntialias)

# Green is reserved, everywhere and in every view, for metabolites that are NOT native to
# the model — compounds a designed pathway introduced. That distinction is the one a user
# acts on, so no compartment may borrow the hue however many compartments a model has.
_RESERVED_HUE_LOW, _RESERVED_HUE_HIGH = 95, 165

# Distinct hue per compartment so compartments are separable at a glance; the
# cytosol keeps the base metabolite blue, others get their own steady hue.
_COMPARTMENT_HUES = {
    "c": None, "cytosol": None,          # base blue
    "e": 28, "extracellular": 28,        # orange
    "p": 200, "periplasm": 200,          # cyan
    "m": 275, "mitochondria": 275,       # purple
    "x": 320, "peroxisome": 320,         # magenta (was green — now reserved)
    "r": 340, "h": 55, "g": 45, "n": 250, "v": 190, "l": 12,
}
_COMP_CACHE: Dict[str, str] = {}


def _avoid_reserved_hue(hue: int) -> int:
    """Push a hue out of the green band kept for non-native metabolites."""
    if _RESERVED_HUE_LOW <= hue <= _RESERVED_HUE_HIGH:
        span = _RESERVED_HUE_HIGH - _RESERVED_HUE_LOW
        # Rotate past the band rather than clamping to its edge, so two compartments
        # that both land inside it do not collapse onto the same colour.
        return (_RESERVED_HUE_HIGH + (hue - _RESERVED_HUE_LOW) * (360 - span) // span) % 360
    return hue


def _compartment_shade(comp: str) -> str:
    """A metabolite colour shaded by compartment (base blue for the cytosol)."""
    comp = (comp or "").lower()
    if comp in _COMP_CACHE:
        return _COMP_CACHE[comp]
    base = QColor(style.NODE_METABOLITE)
    hue = _COMPARTMENT_HUES.get(comp, None)
    if hue is None and comp not in ("c", "cytosol", ""):
        # Unknown compartment: derive a stable hue from its name.
        hue = (sum(ord(ch) for ch in comp) * 47) % 360
    if hue is None:
        col = base
    else:
        col = QColor.fromHsv(int(_avoid_reserved_hue(int(hue))), 150, 210)
    _COMP_CACHE[comp] = col.name()
    return _COMP_CACHE[comp]


def compartment_legend(model) -> List[Tuple[str, str, str]]:
    """(compartment id, readable name, colour) for every compartment a model uses.

    Model-dependent by construction: a three-compartment model gets three entries. Only
    compartments that actually carry metabolites are listed, because a legend entry for
    something not on screen is worse than no entry at all.
    """
    if model is None:
        return []
    used = {(met.compartment or "").strip() for met in model.metabolites}
    used.discard("")
    names = dict(getattr(model, "compartments", {}) or {})
    out = []
    for comp in sorted(used):
        label = (names.get(comp) or "").strip()
        out.append((comp, label or comp, _compartment_shade(comp)))
    return out


class _NodeItem(QGraphicsItem):
    """A custom-painted, draggable metabolite/reaction node.

    Reaction nodes vary their icon by ``rxn_type`` (reversible / irreversible /
    transport / exchange). Connected edge lines are kept attached as the node is
    dragged.
    """

    def __init__(self, kind: str, label: str, color: str, size: float, rxn_type: str = ""):
        super().__init__()
        self.kind = kind
        self.rxn_type = rxn_type
        self._color = QColor(color)
        self._size = size
        self._hover = False
        self._highlight = False
        self.edges: List[Tuple["_EdgeItem", int]] = []
        self.setData(_DATA_LABEL, label)
        self.setData(_DATA_KIND, kind)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.OpenHandCursor)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemSendsScenePositionChanges, True)
        self.setZValue(1)

    def boundingRect(self) -> QRectF:  # noqa: N802
        s = self._size + 6
        return QRectF(-s, -s, 2 * s, 2 * s)

    def set_highlighted(self, on: bool) -> None:
        if on != self._highlight:
            self._highlight = on
            self.update()

    # keep connected edges attached while dragging
    def itemChange(self, change, value):  # noqa: N802
        if change == QGraphicsItem.ItemScenePositionHasChanged:
            p = self.scenePos()
            for edge, end in self.edges:
                edge.set_endpoint(end, p.x(), p.y())
        return super().itemChange(change, value)

    def hoverEnterEvent(self, event):  # noqa: N802
        self._hover = True
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):  # noqa: N802
        self._hover = False
        self.update()
        super().hoverLeaveEvent(event)

    def paint(self, painter: QPainter, option, widget=None):  # noqa: N802
        painter.setRenderHint(QPainter.Antialiasing, True)
        s = self._size

        if self._hover or self._highlight:
            halo = QColor(style.NODE_HIGHLIGHT if self._highlight else self._color)
            halo.setAlpha(70 if self._hover and not self._highlight else 110)
            painter.setPen(Qt.NoPen)
            painter.setBrush(QBrush(halo))
            painter.drawEllipse(QPointF(0, 0), s + 5, s + 5)

        outline = QPen(QColor("#2b2f33"), max(0.8, s * 0.06))
        if self.kind == "metabolite":
            self._paint_molecule(painter, s, outline)
        elif self.rxn_type == "exchange":
            self._paint_badge(painter, s, outline, "E")
        elif self.rxn_type == "transport":
            self._paint_badge(painter, s, outline, "T", membrane=True)
        else:
            self._paint_enzyme(painter, s, outline, reversible=self.rxn_type == "reversible")

    def _paint_molecule(self, painter: QPainter, s: float, outline: QPen) -> None:
        center = self._color
        satellite = QColor(center).lighter(135)
        painter.setPen(QPen(QColor("#9aa6c2"), max(0.8, s * 0.12)))
        for ang in (90, 210, 330):
            r = math.radians(ang)
            painter.drawLine(QPointF(0, 0), QPointF(math.cos(r) * s * 0.85, -math.sin(r) * s * 0.85))
        painter.setPen(outline)
        painter.setBrush(QBrush(satellite))
        for ang in (90, 210, 330):
            r = math.radians(ang)
            c = QPointF(math.cos(r) * s * 0.85, -math.sin(r) * s * 0.85)
            painter.drawEllipse(c, s * 0.34, s * 0.34)
        painter.setBrush(QBrush(center))
        painter.drawEllipse(QPointF(0, 0), s * 0.6, s * 0.6)

    def _paint_enzyme(self, painter: QPainter, s: float, outline: QPen, reversible: bool) -> None:
        hexagon = QPolygonF()
        for i in range(6):
            ang = math.radians(60 * i + 30)
            hexagon.append(QPointF(math.cos(ang) * s, math.sin(ang) * s))
        painter.setPen(outline)
        painter.setBrush(QBrush(self._color))
        painter.drawPolygon(hexagon)

        pen = QPen(QColor("white"), max(1.0, s * 0.16))
        painter.setPen(pen)
        painter.drawLine(QPointF(-s * 0.5, 0), QPointF(s * 0.5, 0))
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(QColor("white")))
        # forward arrowhead
        painter.drawPolygon(QPolygonF([
            QPointF(s * 0.2, -s * 0.3), QPointF(s * 0.55, 0), QPointF(s * 0.2, s * 0.3)]))
        if reversible:
            # backward arrowhead -> double-headed (↔)
            painter.drawPolygon(QPolygonF([
                QPointF(-s * 0.2, -s * 0.3), QPointF(-s * 0.55, 0), QPointF(-s * 0.2, s * 0.3)]))

    def _paint_badge(self, painter: QPainter, s: float, outline: QPen, letter: str,
                     membrane: bool = False) -> None:
        rect = QRectF(-s, -s * 0.8, 2 * s, 1.6 * s)
        painter.setPen(outline)
        painter.setBrush(QBrush(self._color))
        painter.drawRoundedRect(rect, s * 0.35, s * 0.35)
        if membrane:
            # two membrane rails to evoke a transporter crossing a bilayer
            painter.setPen(QPen(QColor(255, 255, 255, 150), max(0.8, s * 0.1)))
            painter.drawLine(QPointF(-s, -s * 0.35), QPointF(s, -s * 0.35))
            painter.drawLine(QPointF(-s, s * 0.35), QPointF(s, s * 0.35))
        font = QFont("Segoe UI", 1)
        font.setBold(True)
        font.setPointSizeF(max(4.0, s * 1.05))
        painter.setFont(font)
        painter.setPen(QPen(QColor("white")))
        painter.drawText(rect, Qt.AlignCenter, letter)


class _EdgeItem(QGraphicsItem):
    """A reaction↔metabolite connector that paints the line AND (optionally) an
    arrowhead in one item, so both follow the nodes when dragged — no stale ghosts.

    ``arrow`` is None (plain), or 'p1'/'p2' naming the endpoint the arrow points at
    (irreversible flow); ``dashed`` marks transport reactions."""

    def __init__(self, x1, y1, x2, y2, color, width, dashed=False, arrow=None,
                 arrow_size=9.0):
        super().__init__()
        self._p1 = QPointF(x1, y1)
        self._p2 = QPointF(x2, y2)
        self._color = QColor(color)
        self._width = float(width)
        self._dashed = dashed
        self._arrow = arrow
        self._arrow_size = float(arrow_size)
        self.setZValue(-1)

    def set_endpoint(self, index: int, x: float, y: float) -> None:
        self.prepareGeometryChange()
        if index == 0:
            self._p1 = QPointF(x, y)
        else:
            self._p2 = QPointF(x, y)
        self.update()

    def boundingRect(self) -> QRectF:  # noqa: N802
        extra = self._width + self._arrow_size + 8
        return QRectF(self._p1, self._p2).normalized().adjusted(-extra, -extra, extra, extra)

    def paint(self, painter: QPainter, option, widget=None):  # noqa: N802
        painter.setRenderHint(QPainter.Antialiasing, True)
        pen = QPen(self._color, self._width)
        if self._dashed:
            pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        painter.drawLine(self._p1, self._p2)
        if not self._arrow:
            return
        tail = self._p2 if self._arrow == "p1" else self._p1
        head = self._p1 if self._arrow == "p1" else self._p2
        dx, dy = head.x() - tail.x(), head.y() - tail.y()
        length = math.hypot(dx, dy)
        if length < 12:
            return
        ux, uy = dx / length, dy / length
        tip = QPointF(tail.x() + ux * length * 0.6, tail.y() + uy * length * 0.6)
        # Scale the arrowhead with the node size so it stays clearly visible next
        # to the metabolite/reaction icons, especially on small networks (#B8).
        b = max(9.0, self._arrow_size)
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(self._color))
        painter.drawPolygon(QPolygonF([
            tip,
            QPointF(tip.x() - ux * b - uy * b * 0.55, tip.y() - uy * b + ux * b * 0.55),
            QPointF(tip.x() - ux * b + uy * b * 0.55, tip.y() - uy * b - ux * b * 0.55)]))


class _DetailsPanel(QWidget):
    """Side panel showing the selected node's properties (+ structure / flux)."""

    def __init__(self):
        super().__init__()
        self.setMinimumWidth(220)
        self._fetcher: Optional[StructureFetcher] = None
        self._current_met: Optional[str] = None

        self._title = QLabel("Select a node")
        self._title.setStyleSheet("font-weight:700; font-size:13px; padding:4px 0;")
        self._title.setWordWrap(True)

        self._form_host = QWidget()
        self._form = QFormLayout(self._form_host)
        self._form.setLabelAlignment(Qt.AlignLeft)

        self._structure = QLabel()
        self._structure.setAlignment(Qt.AlignCenter)
        self._structure.setMinimumHeight(10)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.addWidget(self._title)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setStyleSheet("color:#E0E3E8;")
        lay.addWidget(line)
        lay.addWidget(self._form_host)
        lay.addWidget(self._structure)
        lay.addStretch(1)
        self.show_placeholder()

    def _clear_form(self) -> None:
        while self._form.count():
            item = self._form.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._structure.clear()

    def _row(self, key: str, value: str) -> None:
        v = QLabel(str(value))
        v.setWordWrap(True)
        v.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self._form.addRow(f"{key}:", v)

    def show_placeholder(self) -> None:
        self._title.setText("Select a node")
        self._clear_form()
        self._row("Tip", "Click a reaction or metabolite in the map.")

    def show_reaction(self, rxn: cobra.Reaction, flux=None) -> None:
        self._current_met = None
        self._title.setText(f"Reaction · {rxn.id}")
        self._clear_form()
        ec = reaction_ec_numbers(rxn)
        direction = "reversible (↔)" if rxn.reversibility else "irreversible (→)"
        if rxn.boundary:
            direction = "exchange / boundary"
        self._row("Name", network_graph.display_reaction_name(rxn) or "—")
        self._row("Equation", network_graph.reaction_equation(rxn))
        self._row("Direction", direction)
        self._row("Bounds", f"[{rxn.lower_bound:g}, {rxn.upper_bound:g}]")
        self._row("EC number", ", ".join(ec) if ec else "—")
        self._row("Gene rule", rxn.gene_reaction_rule or "—")
        if flux is not None:
            self._row("Current flux", f"{flux:.6g}")

    def show_metabolite(self, met: cobra.Metabolite) -> None:
        self._title.setText(f"Metabolite · {met.id}")
        self._clear_form()
        self._row("Name", network_graph.short_metabolite_name(met.id, met.name or "") or "—")
        self._row("Formula", met.formula or "—")
        self._row("Charge", str(met.charge) if met.charge is not None else "—")
        self._row("Compartment", met.compartment or "—")
        for key in ("kegg.compound", "bigg.metabolite", "chebi"):
            val = met.annotation.get(key) if isinstance(met.annotation, dict) else None
            if val:
                self._row(key.split(".")[0].upper(),
                          ", ".join(val) if isinstance(val, (list, tuple)) else str(val))
        self._structure.setText("Loading structure…")
        self._fetch_structure(met)

    def _fetch_structure(self, met: cobra.Metabolite) -> None:
        self._current_met = met.id
        # Prefer drawing from the metabolite's own SMILES/InChI with RDKit (no
        # network); only fall back to a PubChem lookup when those are unavailable.
        if self._fetcher is not None and self._fetcher.isRunning():
            self._fetcher.fetched.disconnect()
        self._fetcher = StructureFetcher.for_metabolite(met, size=300)
        self._fetcher.fetched.connect(self._on_structure)
        self._fetcher.start()

    def _on_structure(self, met_id: str, data: bytes) -> None:
        if met_id != self._current_met:
            return
        if not data:
            self._structure.setText("(2-D structure unavailable)")
            return
        pix = QPixmap()
        if pix.loadFromData(data):
            self._structure.setPixmap(pix.scaledToWidth(220, Qt.SmoothTransformation))
        else:
            self._structure.setText("(2-D structure unavailable)")


class _GraphView(QGraphicsView):
    """Pan/zoom view; drag nodes to reposition; emits click and context signals."""

    node_clicked = Signal(str)
    node_context = Signal(str, str, object)  # label, kind, global QPoint

    def __init__(self):
        super().__init__()
        self.setRenderHint(QPainter.Antialiasing)
        self.setRenderHint(QPainter.TextAntialiasing)
        self.setDragMode(QGraphicsView.NoDrag)
        self.setScene(QGraphicsScene(self))
        self.setBackgroundBrush(QBrush(QColor("#FBFCFD")))
        self.setCursor(Qt.ArrowCursor)
        self.viewport().setCursor(Qt.ArrowCursor)
        self._panning = False
        self._moved = False
        self._last_pan = QPointF()
        self._press_node: Optional[QGraphicsItem] = None
        self._press_pos = QPointF()

    def wheelEvent(self, event):  # noqa: N802
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.scale(factor, factor)

    def _node_item_at(self, pos) -> Optional[QGraphicsItem]:
        item = self.itemAt(pos.toPoint() if hasattr(pos, "toPoint") else pos)
        while item is not None:
            if item.data(_DATA_LABEL):
                return item if isinstance(item, _NodeItem) else item.parentItem() or item
            item = item.parentItem()
        return None

    def mousePressEvent(self, event):  # noqa: N802
        self._press_node = None
        if event.button() == Qt.LeftButton:
            node = self._node_item_at(event.position())
            if isinstance(node, _NodeItem):
                self._press_node = node
                self._press_pos = event.position()
                self._moved = False
                super().mousePressEvent(event)  # let the item start dragging
                return
            self._panning = True
            self._moved = False
            self._last_pan = event.position()
            self.viewport().setCursor(Qt.ClosedHandCursor)
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):  # noqa: N802
        if self._panning:
            delta = event.position() - self._last_pan
            if delta.manhattanLength() > 2:
                self._moved = True
            self._last_pan = event.position()
            self.horizontalScrollBar().setValue(
                int(self.horizontalScrollBar().value() - delta.x()))
            self.verticalScrollBar().setValue(
                int(self.verticalScrollBar().value() - delta.y()))
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):  # noqa: N802
        self.viewport().setCursor(Qt.ArrowCursor)
        if event.button() == Qt.LeftButton and self._press_node is not None:
            moved = (event.position() - self._press_pos).manhattanLength() > 3
            node = self._press_node
            self._press_node = None
            super().mouseReleaseEvent(event)
            if not moved:
                label = node.data(_DATA_LABEL)
                if label:
                    self.node_clicked.emit(label)
            return
        if event.button() == Qt.LeftButton and self._panning:
            self._panning = False
        super().mouseReleaseEvent(event)

    def contextMenuEvent(self, event):  # noqa: N802
        label, kind = self._node_at(event.pos())
        if label:
            self.node_context.emit(label, kind, event.globalPos())

    def _node_at(self, pos) -> Tuple[Optional[str], Optional[str]]:
        item = self.itemAt(pos.toPoint() if hasattr(pos, "toPoint") else pos)
        while item is not None:
            label = item.data(_DATA_LABEL)
            if label:
                return label, item.data(_DATA_KIND)
            item = item.parentItem()
        return None, None


class NetworkView(QWidget):
    node_selected = Signal(str)
    node_context_requested = Signal(str, str, object)  # label, kind, global pos

    LAYOUTS = [("Radial", "radial"), ("Layered (flow)", "layered"), ("Force-directed", "force")]

    def __init__(self):
        super().__init__()
        self._model: Optional[cobra.Model] = None
        self._fluxes: Optional[dict] = None
        self._flux_values: Dict[str, str] = {}   # reaction_id -> flux label text
        self._added_metabolites: set = set()      # metabolites added to the model (green)
        self._categories: List = []
        self._category_colors: Dict[str, str] = {}
        self._highlight: Optional[str] = None
        self._node_items: Dict[str, _NodeItem] = {}

        self._subsystems: Dict[str, list] = {}     # {subsystem: [reaction ids]}
        self._selected_scope: set = set()          # chosen category / subsystem names

        # Group focus: Whole model / Category / Subsystem, with a multi-select of one
        # or several categories or subsystems (same as the Escher Visualizer).
        self.focus_kind = QComboBox()
        self.focus_kind.addItems(["Whole model", "Category", "Subsystem"])
        self.focus_kind.setToolTip("Draw the whole model, or focus on one/several "
                                   "categories or subsystems for a clean view.")
        self.focus_kind.currentIndexChanged.connect(self._on_focus_kind_changed)
        self.scope_btn = QPushButton("Select…")
        self.scope_btn.clicked.connect(self._open_scope_selector)
        self.scope_btn.setVisible(False)

        self.focus_combo = QComboBox()
        self.focus_combo.setEditable(True)
        self.focus_combo.setInsertPolicy(QComboBox.NoInsert)
        # Bounded width: long reaction/metabolite names must not dictate the view's
        # (and therefore the whole window's) minimum width. The drop-down popup is
        # still as wide as its contents.
        self.focus_combo.setMinimumWidth(150)
        self.focus_combo.setMaximumWidth(260)
        self.focus_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.focus_combo.setMinimumContentsLength(12)
        self.focus_combo.setToolTip(
            "Center the map on a reaction/metabolite. Type to search, or click the arrow "
            "to pick from the list. (Use the group focus on the left for whole categories "
            "or subsystems.)")
        self.focus_combo.activated.connect(self.redraw)

        self.radius_spin = QSpinBox()
        self.radius_spin.setRange(0, 6)
        self.radius_spin.setValue(1)
        self.radius_spin.setToolTip(
            "How many reaction steps to expand around the focus. For a category, 0 shows just "
            "the category; 1+ adds reactions that many steps away across the whole network.")

        self.layout_combo = QComboBox()
        for label, key in self.LAYOUTS:
            self.layout_combo.addItem(label, userData=key)
        self.layout_combo.setToolTip(
            "Radial: rings around the focus. Layered: left-to-right pathway flow. "
            "Force-directed: organic overview.")
        self.layout_combo.currentIndexChanged.connect(self.redraw)

        self.currency_check = QCheckBox("Hide currency metabolites")
        self.currency_check.setChecked(True)
        self.currency_check.setToolTip(
            "Hide ubiquitous metabolites (ATP, H2O, H+, NAD(P)H …) that clutter the map.")

        draw_btn = QPushButton("Draw")
        draw_btn.clicked.connect(self.redraw)
        export_btn = QPushButton("Export PNG…")
        export_btn.setToolTip("Save the current network map as a high-resolution PNG image.")
        export_btn.clicked.connect(self.export_png)
        self.details_btn = QPushButton("Details ▸")
        self.details_btn.setCheckable(True)
        self.details_btn.setChecked(True)
        self.details_btn.setToolTip("Show/hide the details panel for the selected node.")
        self.details_btn.toggled.connect(self._toggle_details)

        # A wrapping toolbar: it reflows onto extra rows when narrow instead of
        # cropping, and — crucially — its minimum width is just the widest single
        # control, so this view can never force the main window wider than the screen.
        from ..widgets.flow_layout import FlowLayout
        toolbar = FlowLayout()
        toolbar.addWidget(QLabel("Focus:"))
        toolbar.addWidget(self.focus_kind)
        toolbar.addWidget(self.scope_btn)
        toolbar.addWidget(QLabel("Center on:"))
        toolbar.addWidget(self.focus_combo)
        toolbar.addWidget(QLabel("Steps:"))
        toolbar.addWidget(self.radius_spin)
        toolbar.addWidget(QLabel("Layout:"))
        toolbar.addWidget(self.layout_combo)
        toolbar.addWidget(self.currency_check)
        toolbar.addWidget(draw_btn)
        toolbar.addWidget(export_btn)
        toolbar.addWidget(self.details_btn)

        self.view = _GraphView()
        self.view.node_clicked.connect(self._on_node_clicked)
        self.view.node_context.connect(self.node_context_requested)

        self.status = QLabel("Open a model to draw its network.")
        self.status.setStyleSheet(f"color: {style.TEXT_MUTED}; padding: 2px;")

        graph_side = QWidget()
        graph_layout = QVBoxLayout(graph_side)
        graph_layout.setContentsMargins(0, 0, 0, 0)
        graph_layout.addWidget(self.view, 1)
        graph_layout.addWidget(self.status)

        self.details = _DetailsPanel()
        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.addWidget(graph_side)
        self.splitter.addWidget(self.details)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([900, 280])

        # Which colour means what, for the model actually open. Built from the model
        # rather than fixed, so a three-compartment model shows three entries.
        self.legend = QLabel()
        self.legend.setTextFormat(Qt.RichText)
        self.legend.setWordWrap(True)
        self.legend.setStyleSheet("padding: 2px 4px;")
        self._update_legend()

        layout = QVBoxLayout(self)
        layout.addLayout(toolbar)
        layout.addWidget(self.legend)
        layout.addWidget(self.splitter, 1)

    def _update_legend(self) -> None:
        """Restate what the colours mean for the model currently open."""
        def swatch(colour: str, text: str) -> str:
            return (f"<span style='color:{colour};font-size:13px;'>&#9679;</span>&nbsp;"
                    f"{text}")

        parts = [swatch(style.NODE_ADDED_METABOLITE, "<b>not native to the model</b>")]
        for comp, name, colour in compartment_legend(self._model):
            label = name if name.lower() != comp.lower() else f"compartment “{comp}”"
            parts.append(swatch(colour, f"{label} ({comp})"))
        parts.append(swatch(style.NODE_REACTION, "reaction"))

        if len(parts) <= 2:
            self.legend.setText(
                f"<span style='color:{style.TEXT_MUTED};'>"
                f"{swatch(style.NODE_ADDED_METABOLITE, 'not native to the model')}"
                f"</span>")
            return
        self.legend.setText(
            f"<span style='color:{style.TEXT_MUTED};'>"
            + " &nbsp;&nbsp;".join(parts) + "</span>")

    # ----- public API --------------------------------------------------
    def set_model(self, model: cobra.Model) -> None:
        self._model = model
        self._fluxes = None
        self._update_legend()
        self._highlight = None
        from ..widgets.scope_dialog import model_subsystems
        self._subsystems = model_subsystems(model)
        self._selected_scope = set()
        self._update_scope_button()
        self.details.show_placeholder()
        self._populate_focus()
        self._redraw_or_defer()

    # ----- group focus (multiple categories / subsystems) -------------------
    def _on_focus_kind_changed(self) -> None:
        kind = self.focus_kind.currentText()
        self._selected_scope = set()
        self.scope_btn.setVisible(kind in ("Category", "Subsystem"))
        self._update_scope_button()
        self.redraw()

    def _scope_items(self) -> dict:
        kind = self.focus_kind.currentText()
        if kind == "Category":
            return {c.name: sorted(c.reaction_ids) for c in self._categories}
        if kind == "Subsystem":
            return self._subsystems
        return {}

    def _update_scope_button(self) -> None:
        kind = self.focus_kind.currentText()
        if kind == "Whole model":
            self.scope_btn.setVisible(False)
            return
        self.scope_btn.setVisible(True)
        noun = "categories" if kind == "Category" else "subsystems"
        n = len(self._selected_scope)
        self.scope_btn.setText(f"Select {noun}…" if n == 0
                               else f"{n} selected")

    def _open_scope_selector(self) -> None:
        from PySide6.QtWidgets import QMessageBox

        from ..widgets.scope_dialog import choose_scope
        kind = self.focus_kind.currentText()
        items = self._scope_items()
        if not items:
            QMessageBox.information(
                self, f"No {kind.lower()}s",
                f"This model has no {kind.lower()}s defined."
                + (" Create categories in the Categories panel first."
                   if kind == "Category" else ""))
            return
        chosen = choose_scope(self, kind, items, self._selected_scope)
        if chosen is None:
            return
        self._selected_scope = chosen
        self._update_scope_button()
        self.redraw()

    def _group_focus_ids(self):
        """Union of reaction ids for the selected categories/subsystems, or None."""
        if self.focus_kind.currentText() == "Whole model" or not self._selected_scope:
            return None
        items = self._scope_items()
        ids: set = set()
        for name in self._selected_scope:
            ids.update(items.get(name, []))
        return ids or None

    def set_categories(self, categories) -> None:
        self._categories = list(categories)
        self._category_colors = {}
        for cat in self._categories:
            for rid in cat.reaction_ids:
                self._category_colors[rid] = cat.color
        if self.focus_kind.currentText() == "Category":
            self._selected_scope &= {c.name for c in self._categories}
            self._update_scope_button()
        self._populate_focus()
        self._redraw_or_defer()

    def set_fluxes(self, fluxes: Optional[dict]) -> None:
        self._fluxes = fluxes
        self._redraw_or_defer()

    def set_flux_values(self, values: Optional[dict]) -> None:
        """Text (e.g. '12.3' or '[0, 5.2]') shown under each reaction node."""
        self._flux_values = dict(values or {})
        self._redraw_or_defer()

    def set_added_metabolites(self, ids) -> None:
        """Metabolite ids not in the original model — drawn green (vs native blue)."""
        self._added_metabolites = set(ids or ())
        self._redraw_or_defer()

    def set_manual_draw(self, on: bool) -> None:
        """When on, new data NEVER auto-renders — the user must click Draw. Used by
        the Strategy Visualizer, where laying out a genome-scale difference map the
        moment a second strategy is saved is slow and rarely what the user wants."""
        self._manual_draw = bool(on)

    def _redraw_or_defer(self) -> None:
        """Render now if the map is on screen; otherwise defer until it is shown.

        Laying out a genome-scale network is expensive, so we don't do it while the
        Network Map tab is hidden (e.g. right after loading a model, when the
        Analysis tab is in front) — the app becomes usable immediately and the map
        draws the first time the user opens the tab. In manual-draw mode nothing is
        rendered until the user asks for it."""
        if getattr(self, "_manual_draw", False):
            self._needs_redraw = False
            return
        if self.isVisible():
            self.redraw()
        else:
            self._needs_redraw = True

    def showEvent(self, event) -> None:  # noqa: N802 - Qt override
        super().showEvent(event)
        if getattr(self, "_manual_draw", False):
            return
        if getattr(self, "_needs_redraw", False):
            self._needs_redraw = False
            self._redraw_with_feedback()

    def set_render_busy(self, on: bool) -> None:
        """Whether the first render shows a modal 'building the map' popup. Off for
        dialog-embedded views (focused sub-networks render fast, and a nested modal
        popup could get stuck behind the dialog's own event loop — #B8)."""
        self._render_busy = on

    def _redraw_with_feedback(self) -> None:
        """Render the map, showing a loading popup for non-trivial models so the
        first open of the Network Map tab doesn't look like a frozen window."""
        big = self._model is not None and len(self._model.reactions) > 60
        if big and getattr(self, "_render_busy", True):
            from ..widgets.busy import run_on_gui_with_busy
            run_on_gui_with_busy(self.window(), "Building the network map…", self.redraw,
                                 title="Network map")
        else:
            self.redraw()

    def focus_on(self, node_id: str, steps: int = 1) -> None:
        if self._model is None:
            return
        idx = self.focus_combo.findData(node_id)
        if idx >= 0:
            self.focus_combo.setCurrentIndex(idx)
        self.radius_spin.setValue(max(1, steps))
        self._highlight = self._node_key(node_id)
        self.redraw()

    def focus_category(self, name: str) -> None:
        idx = self.focus_combo.findData(f"cat::{name}")
        if idx >= 0:
            self.focus_combo.setCurrentIndex(idx)
        self._highlight = None
        self.redraw()

    def highlight_node(self, node_id: str) -> None:
        key = self._node_key(node_id)
        self._highlight = key
        for nid, item in self._node_items.items():
            item.set_highlighted(nid == key)

    def export_png(self) -> None:
        scene = self.view.scene()
        rect = scene.itemsBoundingRect().adjusted(-40, -40, 40, 40)
        if not rect.isValid() or rect.isEmpty():
            return
        from ..widgets.dialog_util import choose_save_path
        path = choose_save_path(self, "Export network map", "network_map.png",
                                "PNG image (*.png)")
        if not path:
            return
        if not path.lower().endswith(".png"):
            path += ".png"
        scale = 2.0
        image = QImage(int(rect.width() * scale), int(rect.height() * scale),
                       QImage.Format_ARGB32)
        image.fill(QColor("#FFFFFF"))
        painter = QPainter(image)
        painter.setRenderHint(QPainter.Antialiasing, True)
        scene.render(painter, target=QRectF(image.rect()), source=rect)
        painter.end()
        image.save(path)
        self.status.setText(f"Saved network map to {path}")

    # ----- focus combo population --------------------------------------
    def _populate_focus(self) -> None:
        if self._model is None:
            return
        prev = self.focus_combo.currentData()
        self.focus_combo.blockSignals(True)
        self.focus_combo.clear()
        self.focus_combo.addItem("Whole model", userData=None)
        for cat in self._categories:
            self.focus_combo.addItem(f"◆ Category: {cat.name}", userData=f"cat::{cat.name}")
        for rxn in list(self._model.reactions)[:5000]:
            self.focus_combo.addItem(f"R: {rxn.id}", userData=rxn.id)
        for met in list(self._model.metabolites)[:5000]:
            self.focus_combo.addItem(f"M: {met.id}", userData=met.id)
        idx = self.focus_combo.findData(prev) if prev is not None else 0
        self.focus_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.focus_combo.blockSignals(False)

    # ----- rendering ---------------------------------------------------
    def redraw(self) -> None:
        scene = self.view.scene()
        scene.clear()
        self._node_items.clear()
        if self._model is None:
            return

        center, subset, seeds = self._resolve_focus()
        try:
            graph = network_graph.build_graph(
                self._model,
                center=center,
                radius=self.radius_spin.value(),
                hide_currency=self.currency_check.isChecked(),
                fluxes=self._fluxes,
                reaction_subset=subset,
                seed_reactions=seeds,
                layout=self.layout_combo.currentData(),
            )
        except KeyError as exc:
            self.status.setText(str(exc))
            return

        self._render(graph)
        flux_note = " · flux overlay on" if self._fluxes else ""
        if seeds is not None:
            scope = f"category +{self.radius_spin.value()} step(s)"
        elif subset is not None:
            scope = "category"
        else:
            scope = "view"
        self.status.setText(
            f"{graph.node_count()} nodes in {scope}{flux_note}. Click to inspect · drag nodes "
            "to arrange · right-click for actions · scroll to zoom · drag canvas to pan.")

    def _resolve_focus(self):
        """Return (center_id, subset_ids, seed_ids) from the current Focus selection."""
        # Group focus (one/several categories or subsystems) takes precedence: Steps 0
        # shows exactly that group; Steps > 0 expands it across the network.
        group = self._group_focus_ids()
        if group is not None:
            if self.radius_spin.value() <= 0:
                return None, set(group), None
            return None, None, set(group)
        data = self.focus_combo.currentData()
        if isinstance(data, str) and data.startswith("cat::"):
            name = data[len("cat::"):]
            ids = set()
            for cat in self._categories:
                if cat.name == name:
                    ids = set(cat.reaction_ids)
                    break
            # Steps 0 -> just the category; >0 -> expand over the whole network.
            if self.radius_spin.value() <= 0:
                return None, ids, None
            return None, None, ids
        if isinstance(data, str):
            return data, None, None
        typed = self.focus_combo.currentText().split(":")[-1].strip()
        if typed and (self._model.reactions.has_id(typed) or self._model.metabolites.has_id(typed)):
            return typed, None, None
        return None, None, None

    def _node_size(self, n_nodes: int) -> float:
        if n_nodes <= 15:
            return 26.0
        if n_nodes <= 40:
            return 18.0
        if n_nodes <= 90:
            return 12.0
        if n_nodes <= 180:
            return 8.0
        return 6.0

    def _render(self, graph: network_graph.NetworkGraph) -> None:
        scene = self.view.scene()
        positions = {n.node_id: QPointF(n.x * _SCALE, n.y * _SCALE) for n in graph.nodes}
        size = self._node_size(graph.node_count())

        max_flux = 0.0
        if self._fluxes:
            max_flux = max((abs(v) for v in self._fluxes.values()), default=0.0)

        # Nodes first so edges can attach to them for drag-tracking.
        for node in graph.nodes:
            if node.kind == "metabolite":
                if node.label in self._added_metabolites:
                    color = style.NODE_ADDED_METABOLITE      # green: added to the model
                else:
                    color = _compartment_shade(node.data.get("compartment", ""))
            else:
                color = self._category_colors.get(node.label, style.NODE_REACTION)
            item = _NodeItem(node.kind, node.label, color, size,
                             rxn_type=node.data.get("rxn_type", ""))
            item.setPos(positions[node.node_id])
            item.setToolTip(f"{node.label}\n{node.data.get('name', '')}")
            item.set_highlighted(node.node_id == self._highlight)
            scene.addItem(item)
            self._node_items[node.node_id] = item
            label_text = node.data.get("display", node.label)
            flux_text = self._flux_values.get(node.label, "") if node.kind == "reaction" else ""
            # When flux values are shown, bold the name and put the value on its own line.
            self._add_label(item, label_text, size, sublabel=flux_text, bold=bool(flux_text))

        rxn_types = {n.node_id: n.data.get("rxn_type", "")
                     for n in graph.nodes if n.kind == "reaction"}
        for edge in graph.edges:
            p1, p2 = positions.get(edge.source), positions.get(edge.target)
            if p1 is None or p2 is None:
                continue
            # Identify the reaction endpoint + its type; the metabolite is the other.
            if edge.source in rxn_types:
                rxn_node, met_pt, rxn_pt = edge.source, p2, p1
            else:
                rxn_node, met_pt, rxn_pt = edge.target, p1, p2
            rtype = rxn_types.get(rxn_node, "")

            # Scale line thickness with the node icon size so connectors stay
            # proportional to the metabolite/reaction icons — thin lines next to
            # big icons were nearly invisible on small networks (#A2). Flux
            # magnitude then modulates on top of that icon-proportional base.
            base_w = max(1.4, size * 0.30)
            color = QColor(style.EDGE)
            width = base_w
            if edge.flux is not None and max_flux > 0:
                magnitude = abs(edge.flux) / max_flux
                color = QColor(style.FLUX_FORWARD if edge.flux >= 0 else style.FLUX_REVERSE)
                width = base_w * (0.7 + 2.8 * magnitude)
            # Irreversible reactions: arrowhead in the direction of flow.
            arrow = None
            if rtype == "irreversible":
                head_is_met = edge.stoichiometry > 0     # metabolite is a product
                head_pt = met_pt if head_is_met else rxn_pt
                arrow = "p1" if head_pt is p1 else "p2"
            item = _EdgeItem(p1.x(), p1.y(), p2.x(), p2.y(), color, width,
                             dashed=(rtype == "transport"), arrow=arrow,
                             arrow_size=max(9.0, size * 1.6))
            scene.addItem(item)
            src_item = self._node_items.get(edge.source)
            tgt_item = self._node_items.get(edge.target)
            if src_item is not None:
                src_item.edges.append((item, 0))
            if tgt_item is not None:
                tgt_item.edges.append((item, 1))

        rect = scene.itemsBoundingRect()
        if rect.isValid():
            self.view.setSceneRect(rect.adjusted(-60, -60, 60, 60))
            self.view.fitInView(self.view.sceneRect(), Qt.KeepAspectRatio)

    def _add_label(self, node_item, text: str, size: float, sublabel: str = "",
                   bold: bool = False) -> None:
        chip = QGraphicsRectItem(node_item)
        chip.setBrush(QBrush(QColor(255, 255, 255, 210)))
        chip.setPen(QPen(Qt.NoPen))
        chip.setZValue(2)

        font = QFont(_LABEL_FONT)
        font.setPointSizeF(max(7.0, size * 0.85))
        font.setBold(bold)                     # emphasise names when flux values shown
        label = QGraphicsSimpleTextItem(text, chip)
        label.setFont(font)
        label.setBrush(QBrush(QColor(style.TEXT)))
        label.setData(_DATA_LABEL, node_item.data(_DATA_LABEL))
        label.setData(_DATA_KIND, node_item.data(_DATA_KIND))

        pad = size * 0.3
        label.setPos(pad, pad * 0.4)
        br = label.boundingRect()
        total_w, total_h = br.width(), br.height()
        if sublabel:
            # A second, non-bold line for the flux value, so it reads apart from the name.
            sub_font = QFont(_LABEL_FONT)
            sub_font.setPointSizeF(max(6.5, size * 0.78))
            sub = QGraphicsSimpleTextItem(sublabel, chip)
            sub.setFont(sub_font)
            sub.setBrush(QBrush(QColor(style.FLUX_FORWARD)))
            sub.setData(_DATA_LABEL, node_item.data(_DATA_LABEL))
            sub.setData(_DATA_KIND, node_item.data(_DATA_KIND))
            sub.setPos(pad, pad * 0.4 + br.height())
            sbr = sub.boundingRect()
            total_w = max(total_w, sbr.width())
            total_h += sbr.height()
        chip.setRect(0, 0, total_w + 2 * pad, total_h + pad * 0.8)
        chip.setPos(size + 2, -size - 2)

    # ----- helpers -----------------------------------------------------
    def _toggle_details(self, on: bool) -> None:
        self.details_btn.setText("Details ▸" if on else "Details ◂")
        self.details.setVisible(on)

    def _node_key(self, node_id: str) -> str:
        if self._model is None:
            return node_id
        if self._model.reactions.has_id(node_id):
            return f"R:{node_id}"
        if self._model.metabolites.has_id(node_id):
            return f"M:{node_id}"
        return node_id

    def _on_node_clicked(self, label: str) -> None:
        self.highlight_node(label)
        self._update_details(label)
        self.node_selected.emit(label)

    def _update_details(self, label: str) -> None:
        if self._model is None:
            return
        if self._model.reactions.has_id(label):
            rxn = self._model.reactions.get_by_id(label)
            flux = self._fluxes.get(label) if self._fluxes else None
            self.details.show_reaction(rxn, flux)
        elif self._model.metabolites.has_id(label):
            self.details.show_metabolite(self._model.metabolites.get_by_id(label))
