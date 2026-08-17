"""Escher Explorer — interactive metabolic maps powered by the Escher engine.

This is the interactive companion to the (matplotlib/Qt) Strategy Explorer. It
embeds the field-standard **Escher** renderer in a QtWebEngine view and drives it
three ways:

  * **Single strategy** — draw one saved strategy's |flux| on the map (sequential
    colour, arrow width ∝ |flux|);
  * **Difference (A − B)** — Δflux between two strategies on a diverging RdBu scale
    (red = more flux after engineering, blue = less);
  * **FBA navigator (live)** — an Escher-FBA-inspired mode: click any reaction on
    the map, tighten/open/knock-out its bounds, and re-solve FBA instantly to watch
    the flux re-route — all powered by the ToolBox's own cobra backend.

Maps are generated on the fly for *any* model (see :mod:`core.escher_map`), so this
works on the user's own models, not just curated BiGG maps. It never overrides the
existing Network Map / Strategy Explorer — it is an additional tab.
"""

from __future__ import annotations

import json
from typing import Dict, Optional

import cobra
from PySide6.QtCore import QObject, Qt, Signal, Slot
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from ...core import escher_map
from ...core.flux_state import StrategyManager

# QtWebEngine is an optional-at-runtime dependency: import lazily so a build
# without it degrades to a clear message instead of failing to start.
try:
    from PySide6.QtCore import QUrl
    from PySide6.QtWebChannel import QWebChannel
    from PySide6.QtWebEngineCore import QWebEnginePage
    from PySide6.QtWebEngineWidgets import QWebEngineView
    _WEBENGINE_OK = True
    _WEBENGINE_ERR = ""
except Exception as _exc:  # noqa: BLE001
    _WEBENGINE_OK = False
    _WEBENGINE_ERR = str(_exc)

_MODES = ["Single strategy", "Difference (A − B)", "FBA navigator (live)"]


class _Bridge(QObject):
    """JS → Python channel. Two DISTINCT interactions on a map element (#5):
      * *selected* — a plain LEFT-click; used to pick a reaction to edit in the FBA
        navigator. Never opens a popup.
      * *info* — chosen from the RIGHT-click menu ("Show reaction details"); opens the
        details popup and nothing else.
    """
    reaction_selected = Signal(str)
    metabolite_selected = Signal(str)
    reaction_info = Signal(str)
    metabolite_info = Signal(str)
    js_ready = Signal()

    @Slot(str)
    def reactionSelected(self, bigg_id: str) -> None:  # noqa: N802 - JS-facing name
        self.reaction_selected.emit(bigg_id)

    @Slot(str)
    def metaboliteSelected(self, bigg_id: str) -> None:  # noqa: N802 - JS-facing name
        self.metabolite_selected.emit(bigg_id)

    @Slot(str)
    def reactionInfo(self, bigg_id: str) -> None:  # noqa: N802 - JS-facing name
        self.reaction_info.emit(bigg_id)

    @Slot(str)
    def metaboliteInfo(self, bigg_id: str) -> None:  # noqa: N802 - JS-facing name
        self.metabolite_info.emit(bigg_id)

    @Slot()
    def jsReady(self) -> None:  # noqa: N802 - JS-facing name
        self.js_ready.emit()


class EscherExplorer(QWidget):
    """Interactive Escher map tab (see module docstring)."""

    # Emitted when the user clicks a reaction/metabolite on the map, so the host
    # window can open its rich details popup (#4/#6).
    reaction_info_requested = Signal(object)
    metabolite_info_requested = Signal(object)

    def __init__(self):
        super().__init__()
        self._model: Optional[cobra.Model] = None
        self._strategies = StrategyManager()
        self._categories = {}                 # {name: [reaction_ids]}
        self._subsystems: Dict[str, list] = {}   # {subsystem: [reaction_ids]}
        self._selected_scope: set = set()     # chosen category / subsystem names
        self._overrides: Dict[str, tuple] = {}   # FBA navigator: rid -> (lb, ub)
        self._selected_rid: Optional[str] = None
        self._page_ready = False
        self._pending_js: list = []

        if not _WEBENGINE_OK:
            self._build_unavailable()
            return
        self._build_ui()

    # ---- fallback when QtWebEngine is missing --------------------------------
    def _build_unavailable(self) -> None:
        lay = QVBoxLayout(self)
        msg = QLabel(
            "<b>Interactive Escher maps need the QtWebEngine component.</b><br><br>"
            "It isn't available in this build, so the interactive map can't be shown "
            "here. The Network Map and Strategy Visualizer tabs still provide flux maps "
            "and difference maps.<br><br>"
            f"<span style='color:#888'>({_WEBENGINE_ERR})</span>")
        msg.setWordWrap(True)
        msg.setAlignment(Qt.AlignCenter)
        lay.addStretch(1)
        lay.addWidget(msg)
        lay.addStretch(1)

    # ---- normal UI -----------------------------------------------------------
    def _build_ui(self) -> None:
        # top control bar
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(_MODES)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        self.combo_a = QComboBox()
        self.combo_b = QComboBox()
        self.ab_dash = QLabel("-")

        # Focus: Whole model / Category / Subsystem. The latter two reveal a button
        # that opens a multi-select of categories or model subsystems (#U4).
        self.focus_kind = QComboBox()
        self.focus_kind.addItems(["Whole model", "Category", "Subsystem"])
        self.focus_kind.setToolTip("Restrict the map to keep it legible. Focusing on "
                                   "one or several subsystems gives the cleanest layout.")
        self.focus_kind.currentIndexChanged.connect(self._on_focus_kind_changed)
        self.scope_btn = QPushButton("Select…")
        self.scope_btn.clicked.connect(self._open_scope_selector)
        self.scope_btn.setVisible(False)

        self.build_btn = QPushButton("Build / refresh map")
        self.build_btn.clicked.connect(self.rebuild)
        # Show/hide the on-map colour legend (#3).
        self.legend_btn = QPushButton("Hide legend")
        self.legend_btn.setCheckable(True)
        self.legend_btn.setChecked(True)          # legend shown by default
        self.legend_btn.setToolTip("Show or hide the colour legend drawn on the map.")
        self.legend_btn.toggled.connect(self._on_legend_toggled)
        self.top_status = QLabel("")
        self.top_status.setStyleSheet("color:#5f6368;")

        # Wrapping toolbar: never crops, and its minimum width is just the widest
        # single control, so this tab can't force the window wider than the screen.
        from ..widgets.flow_layout import FlowLayout
        bar = FlowLayout()
        bar.addWidget(QLabel("View:"))
        bar.addWidget(self.mode_combo)
        bar.addWidget(self.combo_a)
        bar.addWidget(self.ab_dash)
        bar.addWidget(self.combo_b)
        bar.addWidget(QLabel("Focus:"))
        bar.addWidget(self.focus_kind)
        bar.addWidget(self.scope_btn)
        bar.addWidget(self.build_btn)
        bar.addWidget(self.legend_btn)
        bar.addWidget(self.top_status)

        # left column: the FBA-navigator editor (shown only in live navigator mode)
        self.status = self.top_status   # alias: build messages go to the top bar
        self.top_status.setText("Load a model, then Build the map.")

        self.nav_group = QGroupBox("FBA navigator")
        nav = QFormLayout(self.nav_group)
        self.obj_label = QLabel("—")
        self.sel_label = QLabel("<i>Click a reaction on the map to edit its bounds.</i>")
        self.sel_label.setWordWrap(True)
        self.lb_spin = QDoubleSpinBox()
        self.ub_spin = QDoubleSpinBox()
        for s in (self.lb_spin, self.ub_spin):
            s.setRange(-1e6, 1e6)
            s.setDecimals(3)
        self.apply_btn = QPushButton("Apply bounds && re-solve")
        self.apply_btn.clicked.connect(self._apply_selected)
        self.ko_btn = QPushButton("Knock out (0, 0) && re-solve")
        self.ko_btn.clicked.connect(self._knock_out_selected)
        self.reset_btn = QPushButton("Reset all bounds")
        self.reset_btn.clicked.connect(self._reset_overrides)
        nav.addRow("Objective:", self.obj_label)
        nav.addRow(self.sel_label)
        nav.addRow("Lower bound:", self.lb_spin)
        nav.addRow("Upper bound:", self.ub_spin)
        nav.addRow(self.apply_btn)
        nav.addRow(self.ko_btn)
        nav.addRow(self.reset_btn)
        self._set_nav_enabled(False)

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.addWidget(self.nav_group)
        left.addStretch(1)
        self.left_w = QWidget()
        self.left_w.setLayout(left)
        self.left_w.setMinimumWidth(240)
        self.left_w.setMaximumWidth(340)
        self.left_w.setVisible(False)   # only shown in the live FBA navigator (#U6)

        # right: the web map
        self.view = QWebEngineView()
        self.bridge = _Bridge()
        # Left-click SELECTS (navigator editing); right-click ▸ details opens the popup.
        self.bridge.reaction_selected.connect(self._on_reaction_selected)
        self.bridge.metabolite_selected.connect(self._on_metabolite_selected)
        self.bridge.reaction_info.connect(self._on_reaction_info)
        self.bridge.metabolite_info.connect(self._on_metabolite_info)
        self.channel = QWebChannel()
        self.channel.registerObject("bridge", self.bridge)
        self.view.page().setWebChannel(self.channel)
        self.view.loadFinished.connect(self._on_load_finished)
        # Escher's Export SVG/PNG buttons trigger a browser download; capture it and
        # let the user choose where to save (the toolbar in the map, #3).
        try:
            self.view.page().profile().downloadRequested.connect(self._on_download)
        except Exception:  # noqa: BLE001
            pass

        from ...resources import escher_host_html
        self.view.load(QUrl.fromLocalFile(escher_host_html()))

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.left_w)
        splitter.addWidget(self.view)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, 900])

        # One-line hint bar shown only in the live FBA navigator, so the click-to-edit
        # affordance (otherwise invisible) is discoverable (#3.3).
        self.nav_hint = QLabel(
            "🖱  <b>FBA navigator:</b> click any reaction on the map to edit its flux "
            "bounds; right-click ▸ for details.")
        self.nav_hint.setWordWrap(True)
        self.nav_hint.setTextFormat(Qt.RichText)
        self.nav_hint.setStyleSheet(
            "background:#E8F0FE; color:#1967D2; border:1px solid #C6DAFC; "
            "border-radius:4px; padding:4px 8px;")
        self.nav_hint.setVisible(False)

        outer = QVBoxLayout(self)
        outer.addLayout(bar)
        outer.addWidget(self.nav_hint)
        outer.addWidget(splitter, 1)

        self._on_mode_changed()

    # ---- public API (mirrors StrategyExplorer) -------------------------------
    def set_model(self, model: Optional[cobra.Model]) -> None:
        self._model = model
        self._overrides.clear()
        self._selected_rid = None
        # Index the model's subsystems so they can be a Focus option (#U4).
        self._subsystems = {}
        if model is not None:
            for r in model.reactions:
                sub = (getattr(r, "subsystem", "") or "").strip()
                if sub:
                    self._subsystems.setdefault(sub, []).append(r.id)
        self._selected_scope = set()
        if _WEBENGINE_OK:
            self._set_nav_enabled(False)
            self._update_scope_button()

    def set_strategies(self, manager: StrategyManager) -> None:
        self._strategies = manager
        if _WEBENGINE_OK:
            self._refresh_strategy_combos()

    def set_categories(self, categories: dict) -> None:
        """`categories` maps a display name to a list/iterable of reaction ids."""
        self._categories = dict(categories or {})
        # Drop any selected categories that no longer exist.
        if self.focus_kind.currentText() == "Category":
            self._selected_scope &= set(self._categories)
        if _WEBENGINE_OK:
            self._update_scope_button()

    # ---- strategy combo plumbing --------------------------------------------
    def _refresh_strategy_combos(self) -> None:
        names = self._strategies.names()
        for combo in (self.combo_a, self.combo_b):
            cur = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(names)
            idx = combo.findText(cur)
            if idx >= 0:
                combo.setCurrentIndex(idx)
            combo.blockSignals(False)
        if len(names) >= 2 and self.combo_b.currentIndex() == self.combo_a.currentIndex():
            self.combo_b.setCurrentIndex(1)

    def _on_mode_changed(self) -> None:
        idx = self.mode_combo.currentIndex()
        is_diff = idx == 1
        is_nav = idx == 2
        self.combo_a.setVisible(idx in (0, 1))
        self.combo_b.setVisible(is_diff)
        self.ab_dash.setVisible(is_diff)
        # The navigator editor column is shown ONLY in the live FBA navigator; the
        # other modes give the whole window to the map (#U6). The hover/click
        # affordance on reactions is likewise only enabled in navigator mode (#U5).
        self.nav_group.setVisible(is_nav)
        self.left_w.setVisible(is_nav)
        self.nav_hint.setVisible(is_nav)
        self._run_js(f"gsmSetNav({'true' if is_nav else 'false'})")

    # ---- focus / scope selection (#U4) --------------------------------------
    def _on_focus_kind_changed(self) -> None:
        kind = self.focus_kind.currentText()
        self._selected_scope = set()
        self.scope_btn.setVisible(kind in ("Category", "Subsystem"))
        self._update_scope_button()

    def _scope_items(self) -> dict:
        """The available {name: [reaction_ids]} for the current focus kind."""
        kind = self.focus_kind.currentText()
        if kind == "Category":
            return self._categories
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
                               else f"{n} {noun[:-3] if n == 1 else noun} selected")

    def _open_scope_selector(self) -> None:
        kind = self.focus_kind.currentText()
        items = self._scope_items()
        if not items:
            QMessageBox.information(
                self, f"No {kind.lower()}s",
                f"This model has no {kind.lower()}s defined."
                + (" Create categories in the Categories panel first."
                   if kind == "Category" else ""))
            return
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Choose {kind.lower()}(s) to map")
        lay = QVBoxLayout(dlg)
        lay.addWidget(QLabel(f"Tick one or more {kind.lower()}s to draw. Their reactions are "
                             "combined into a single focused map."))
        lst = QListWidget()
        lst.setSelectionMode(QAbstractItemView.NoSelection)
        for name in sorted(items):
            it = QListWidgetItem(f"{name}  ({len(items[name])})")
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if name in self._selected_scope else Qt.Unchecked)
            it.setData(Qt.UserRole, name)
            lst.addItem(it)
        lay.addWidget(lst, 1)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        lay.addWidget(bb)
        dlg.resize(360, 460)
        if dlg.exec() != QDialog.Accepted:
            return
        self._selected_scope = {lst.item(i).data(Qt.UserRole)
                                for i in range(lst.count())
                                if lst.item(i).checkState() == Qt.Checked}
        self._update_scope_button()
        self.rebuild()

    # ---- map building --------------------------------------------------------
    def _scope_reaction_ids(self):
        if self.focus_kind.currentText() == "Whole model" or not self._selected_scope:
            return None
        items = self._scope_items()
        ids: list = []
        seen = set()
        for name in self._selected_scope:
            for rid in items.get(name, []):
                if rid not in seen:
                    seen.add(rid)
                    ids.append(rid)
        return ids or None

    def rebuild(self) -> None:
        if self._model is None:
            self.status.setText("Load a model first.")
            return
        try:
            reaction_ids = self._scope_reaction_ids()
            emap = escher_map.build_escher_map(
                self._model, reaction_ids,
                name=f"{getattr(self._model, 'id', 'model')} — Escher Visualizer")
        except Exception as exc:  # noqa: BLE001
            self.status.setText(f"Could not build the map: {exc}")
            return

        model_json = None
        try:
            model_json = cobra.io.to_json(self._model)
        except Exception:  # noqa: BLE001 - model is optional to Escher
            model_json = None

        mode = self.mode_combo.currentIndex()
        flux, kind = self._current_flux(mode)
        self._run_js(
            f"gsmBuild({json.dumps(json.dumps(emap))}, "
            f"{json.dumps(model_json)}, {json.dumps(json.dumps(flux) if flux is not None else None)}, "
            f"{json.dumps(kind)})")
        if mode == 2:
            self._overrides.clear()
            self._solve_navigator()

    def _current_flux(self, mode: int):
        """Return (flux_dict_or_None, scale_kind) for the selected mode."""
        if mode == 0:  # single strategy — magnitude (abs handled by Escher 'abs' style)
            s = self._strategies.get(self.combo_a.currentText())
            if s is None:
                return None, "single"
            return escher_map.flux_overlay(dict(s.fluxes)), "single"
        if mode == 1:  # difference A − B — signed, diverging (red = up in A vs B)
            a, b = self.combo_a.currentText(), self.combo_b.currentText()
            if not a or not b or a == b:
                self.status.setText("Pick two different strategies for a difference map.")
                return None, "diff"
            # StrategyManager.difference(x, y) returns flux(y) − flux(x); we want
            # flux(A) − flux(B) so positive means "more flux in A" (matches the legend).
            return escher_map.flux_overlay(self._strategies.difference(b, a)), "diff"
        return None, "single"  # navigator seeds flux via _solve_navigator

    # ---- FBA navigator -------------------------------------------------------
    def _set_nav_enabled(self, on: bool) -> None:
        for w in (self.lb_spin, self.ub_spin, self.apply_btn, self.ko_btn):
            w.setEnabled(on)

    def _on_reaction_selected(self, bigg_id: str) -> None:
        """LEFT-click on a reaction: select it. In the FBA navigator this loads it into
        the flux editor. It must NOT open the details popup — that is the whole reason
        info moved to the right-click menu."""
        if self._model is None or not self._model.reactions.has_id(bigg_id):
            return
        rxn = self._model.reactions.get_by_id(bigg_id)
        if self.mode_combo.currentIndex() == 2:
            self._selected_rid = bigg_id
            lb, ub = self._overrides.get(bigg_id, (rxn.lower_bound, rxn.upper_bound))
            from ...core.network_graph import clean_label, display_reaction_name
            self.sel_label.setText(f"<b>{bigg_id}</b> — {display_reaction_name(rxn)}<br>"
                                   f"<span style='color:#888'>"
                                   f"{clean_label(rxn.build_reaction_string())}</span>")
            self.lb_spin.setValue(float(lb))
            self.ub_spin.setValue(float(ub))
            self._set_nav_enabled(True)

    def _on_metabolite_selected(self, bigg_id: str) -> None:
        """LEFT-click on a metabolite: a select gesture. No popup (info is right-click).
        Metabolites have no navigator action, so this is intentionally a no-op beyond
        keeping the interaction symmetric with reactions."""
        return

    def _on_reaction_info(self, bigg_id: str) -> None:
        """RIGHT-click ▸ Show reaction details: open the rich details popup."""
        if self._model is not None and self._model.reactions.has_id(bigg_id):
            self.reaction_info_requested.emit(self._model.reactions.get_by_id(bigg_id))

    def _on_legend_toggled(self, on: bool) -> None:
        """Show/hide the on-map colour legend (#3)."""
        self.legend_btn.setText("Hide legend" if on else "Show legend")
        self._run_js(f"gsmShowLegend({'true' if on else 'false'})")

    def _on_metabolite_info(self, bigg_id: str) -> None:
        """RIGHT-click ▸ Show metabolite details: open the details popup."""
        if self._model is not None and self._model.metabolites.has_id(bigg_id):
            self.metabolite_info_requested.emit(self._model.metabolites.get_by_id(bigg_id))

    def _apply_selected(self) -> None:
        if not self._selected_rid:
            return
        lb, ub = self.lb_spin.value(), self.ub_spin.value()
        if lb > ub:
            lb, ub = ub, lb
        self._overrides[self._selected_rid] = (lb, ub)
        self._solve_navigator()

    def _knock_out_selected(self) -> None:
        if not self._selected_rid:
            return
        self._overrides[self._selected_rid] = (0.0, 0.0)
        self.lb_spin.setValue(0.0)
        self.ub_spin.setValue(0.0)
        self._solve_navigator()

    def _reset_overrides(self) -> None:
        self._overrides.clear()
        self._solve_navigator()

    def _solve_navigator(self) -> None:
        if self._model is None:
            return
        try:
            with self._model as m:
                for rid, (lb, ub) in self._overrides.items():
                    if m.reactions.has_id(rid):
                        r = m.reactions.get_by_id(rid)
                        r.lower_bound, r.upper_bound = lb, ub
                sol = m.optimize()
            flux = escher_map.flux_overlay(dict(sol.fluxes))
            obj = sol.objective_value
            n = len(self._overrides)
            self.obj_label.setText(
                f"{obj:.4g}  ({sol.status}"
                + (f"; {n} edited bound{'s' if n != 1 else ''})" if n else ")"))
            self._run_js(f"gsmSetFlux({json.dumps(json.dumps(flux))})")
        except Exception as exc:  # noqa: BLE001
            self.obj_label.setText(f"infeasible / error: {exc}")

    # ---- web view plumbing ---------------------------------------------------
    def _on_load_finished(self, ok: bool) -> None:
        self._page_ready = bool(ok)
        if ok:
            for js in self._pending_js:
                self.view.page().runJavaScript(js)
            self._pending_js.clear()

    def _run_js(self, js: str) -> None:
        if self._page_ready:
            self.view.page().runJavaScript(js)
        else:
            self._pending_js.append(js)

    def _on_download(self, item) -> None:
        """Save an Escher SVG/PNG export to a user-chosen file (#3)."""
        import os

        from PySide6.QtWidgets import QFileDialog
        try:
            suggested = item.downloadFileName() or "escher_map.svg"
        except Exception:  # noqa: BLE001
            suggested = "escher_map.svg"
        ext = os.path.splitext(suggested)[1].lstrip(".") or "svg"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export map", suggested, f"{ext.upper()} file (*.{ext});;All files (*)")
        if not path:
            item.cancel()
            return
        try:
            item.setDownloadDirectory(os.path.dirname(path))
            item.setDownloadFileName(os.path.basename(path))
            item.accept()
        except Exception:  # noqa: BLE001
            item.cancel()
