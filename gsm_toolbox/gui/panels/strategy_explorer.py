"""Strategy Explorer — the graphical engine's flagship.

A *strategy* is a named solved flux state (one round of engineering). This panel
lists the saved strategies and renders/compares them several ways:

  * Difference map — Δflux between two strategies on the network (diverging colour:
    red = more flux, blue = less), the single most-requested view;
  * Flux map — a single strategy's |flux| on the network;
  * Multi-strategy heatmap — reactions × strategies, values annotated;
  * Titre waterfall — product flux after each round;
  * Parallel coordinates — flux trajectories across strategies.

Strategies come from the project's StrategyManager (saved with the .gsmtbx file).
"""

from __future__ import annotations

import cobra
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from ...core import network_graph
from ...core.flux_state import StrategyManager
from ..views.network_view import NetworkView
from ..viz import plots
from ..viz.plot_canvas import PlotCanvas

_VIEWS = ["Difference map (A vs B)", "Flux map (single)", "Multi-strategy heatmap",
          "Titre waterfall", "Parallel coordinates"]


class StrategyExplorer(QWidget):
    save_strategy_requested = Signal()
    remove_strategy_requested = Signal(str)
    # Run an analysis without leaving this tab (the main window owns the runners).
    run_analysis_requested = Signal(str)     # "fba" | "pfba"

    def __init__(self):
        super().__init__()
        self._model: cobra.Model = None
        self._strategies = StrategyManager()

        # ---- left: strategy list + controls ----
        self.list = QListWidget()
        self.list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list.itemChanged.connect(lambda _i: self._render())
        from ..widgets.wrap_button import WrapButton
        self.save_btn = WrapButton("Save current flux state as strategy…")
        self.save_btn.setToolTip("Capture the current model + last solved fluxes as a named "
                                 "strategy (a round of engineering) to compare later.")
        self.save_btn.clicked.connect(self.save_strategy_requested)
        self.remove_btn = QPushButton("Remove selected")
        self.remove_btn.clicked.connect(self._remove_selected)

        # Run a flux state right here: a strategy captures the LAST solved fluxes, so
        # needing to leave for the Analysis tab and come back is pure friction.
        self.fba_btn = QPushButton("Run FBA")
        self.fba_btn.setToolTip("Solve the current model with FBA, then save the result "
                                "as a strategy.")
        self.fba_btn.clicked.connect(lambda: self.run_analysis_requested.emit("fba"))
        self.pfba_btn = QPushButton("Run pFBA")
        self.pfba_btn.setToolTip("Solve the current model with parsimonious FBA (pFBA), "
                                 "then save the result as a strategy.")
        self.pfba_btn.clicked.connect(lambda: self.run_analysis_requested.emit("pfba"))
        run_row = QHBoxLayout()
        run_row.setContentsMargins(0, 0, 0, 0)
        run_row.addWidget(self.fba_btn)
        run_row.addWidget(self.pfba_btn)
        run_w = QWidget()
        run_w.setLayout(run_row)

        left = QVBoxLayout()
        left.addWidget(QLabel("<b>Strategies</b> (tick to compare)"))
        left.addWidget(self.list, 1)
        left.addWidget(run_w)
        left.addWidget(self.save_btn)
        left.addWidget(self.remove_btn)
        left_w = QWidget()
        left_w.setLayout(left)
        left_w.setMinimumWidth(240)

        # ---- top controls ----
        self.view_combo = QComboBox()
        self.view_combo.addItems(_VIEWS)
        self.view_combo.currentIndexChanged.connect(self._on_view_changed)
        self.combo_a = QComboBox()
        self.combo_b = QComboBox()
        self.combo_a.currentIndexChanged.connect(self._render)
        self.combo_b.currentIndexChanged.connect(self._render)
        self.max_rxn = QSpinBox()
        self.max_rxn.setRange(5, 200)
        self.max_rxn.setValue(25)
        self.max_rxn.setToolTip("Maximum reactions shown in the heatmap / parallel-coordinates "
                                "(most variable across strategies).")
        self.max_rxn.valueChanged.connect(self._render)
        self.ab_label = QLabel("A vs B:")
        # Reaction of interest for the waterfall. The strategy's stored "target" only
        # ever offered exchange reactions, so a heterologous product (or any internal
        # step you actually care about) could not be plotted. Any reaction goes here.
        # Bounded width + editable: an unbounded combo listing thousands of reaction
        # ids would itself widen the window past the screen.
        self.rxn_label = QLabel("Reaction:")
        self.rxn_combo = QComboBox()
        self.rxn_combo.setEditable(True)
        self.rxn_combo.setInsertPolicy(QComboBox.NoInsert)
        self.rxn_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.rxn_combo.setMinimumContentsLength(14)
        self.rxn_combo.setMinimumWidth(160)
        self.rxn_combo.setMaximumWidth(300)
        self.rxn_combo.setToolTip(
            "Reaction to plot across strategies. Defaults to each strategy's saved "
            "target; pick any reaction (including a heterologous step) to follow it.")
        self.rxn_combo.currentIndexChanged.connect(self._render)
        self.rxn_combo.lineEdit().editingFinished.connect(self._render)
        # A wrapping toolbar (never crops, and its minimum width is just the widest
        # single control — so this panel can't force the window wider than the screen).
        from ..widgets.flow_layout import FlowLayout
        controls = FlowLayout()
        controls.addWidget(QLabel("View:"))
        controls.addWidget(self.view_combo)
        controls.addWidget(self.ab_label)
        controls.addWidget(self.combo_a)
        controls.addWidget(self.combo_b)
        controls.addWidget(self.rxn_label)
        controls.addWidget(self.rxn_combo)
        controls.addWidget(QLabel("Max reactions:"))
        controls.addWidget(self.max_rxn)

        # ---- right: stacked map / plot ----
        self.network = NetworkView()
        self.network.set_render_busy(False)   # embedded view (#B8)
        # Never lay out a genome-scale map just because a strategy was saved — the
        # user clicks the map's "Draw" button when they actually want it (#4).
        self.network.set_manual_draw(True)
        self.plot = PlotCanvas()
        self.stack = QStackedWidget()
        self.stack.addWidget(self.network)   # index 0 — map views
        self.stack.addWidget(self.plot)       # index 1 — matplotlib views
        right = QVBoxLayout()
        right.addLayout(controls)
        right.addWidget(self.stack, 1)
        right_w = QWidget()
        right_w.setLayout(right)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_w)
        splitter.addWidget(right_w)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([260, 760])
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(splitter)

        self._placeholder()

    # -- public API -----------------------------------------------------
    def set_model(self, model: cobra.Model) -> None:
        self._model = model
        self.network.set_model(model)
        self._populate_reactions()
        self._render()

    def _populate_reactions(self) -> None:
        """Offer EVERY reaction, exchanges first (the usual product route), so a
        heterologous product added by Pathway Design can be followed too."""
        self.rxn_combo.blockSignals(True)
        cur = self.rxn_combo.currentText()
        self.rxn_combo.clear()
        self.rxn_combo.addItem("(use each strategy's saved target)", "")
        if self._model is not None:
            # cobra's `.exchanges` guesses the external compartment from names and
            # boundary reactions, and RAISES when it cannot — so a model with unusual
            # compartments must not take this panel down with it. Fall back to listing
            # the reactions plainly; ordering is a convenience, not a requirement.
            try:
                ex = sorted(r.id for r in self._model.exchanges)
            except Exception:  # noqa: BLE001
                ex = sorted(r.id for r in self._model.reactions if r.boundary)
            rest = sorted(r.id for r in self._model.reactions
                          if r.id not in set(ex))
            for rid in ex + rest:
                self.rxn_combo.addItem(rid, rid)
            from PySide6.QtWidgets import QCompleter
            comp = QCompleter([self.rxn_combo.itemText(i)
                               for i in range(self.rxn_combo.count())], self.rxn_combo)
            comp.setCaseSensitivity(Qt.CaseInsensitive)
            comp.setFilterMode(Qt.MatchContains)     # type ANY fragment of an id
            self.rxn_combo.setCompleter(comp)
        idx = self.rxn_combo.findText(cur)
        self.rxn_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.rxn_combo.blockSignals(False)

    def _selected_reaction(self) -> str:
        """The reaction to follow, or "" to fall back to each strategy's own target."""
        data = self.rxn_combo.currentData()
        if data:
            return str(data)
        text = self.rxn_combo.currentText().strip()
        if not text or text.startswith("("):
            return ""
        if self._model is not None and self._model.reactions.has_id(text):
            return text
        return ""

    def set_strategies(self, manager: StrategyManager) -> None:
        self._strategies = manager
        self.refresh()

    def refresh(self) -> None:
        names = self._strategies.names()
        self.list.blockSignals(True)
        self.list.clear()
        for n in names:
            it = QListWidgetItem(n)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.list.addItem(it)
        self.list.blockSignals(True)
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
        self.list.blockSignals(False)
        self._render()

    # -- internals ------------------------------------------------------
    def _checked(self) -> list:
        return [self.list.item(i).text() for i in range(self.list.count())
                if self.list.item(i).checkState() == Qt.Checked]

    def _remove_selected(self) -> None:
        for it in self.list.selectedItems():
            self.remove_strategy_requested.emit(it.text())

    def _on_view_changed(self) -> None:
        is_diff = self.view_combo.currentIndex() == 0
        is_single = self.view_combo.currentIndex() == 1
        self.ab_label.setVisible(is_diff or is_single)
        self.combo_a.setVisible(is_diff or is_single)
        self.combo_b.setVisible(is_diff)
        is_waterfall = self.view_combo.currentIndex() == 3
        self.rxn_label.setVisible(is_waterfall)
        self.rxn_combo.setVisible(is_waterfall)
        self.max_rxn.setVisible(self.view_combo.currentIndex() in (2, 4))
        self._render()

    def _label_map(self, ids) -> dict:
        out = {}
        for rid in ids:
            if self._model is not None and self._model.reactions.has_id(rid):
                r = self._model.reactions.get_by_id(rid)
                out[rid] = network_graph.short_metabolite_name(rid, r.name or "") or rid
            else:
                out[rid] = rid
        return out

    def _placeholder(self) -> None:
        self.plot.render(plots._empty,
                         "Save a flux state as a strategy (after running FBA/pFBA), then "
                         "compare strategies here.", title="strategy_explorer")
        self.stack.setCurrentWidget(self.plot)

    def _top_reactions(self, names, limit) -> list:
        """Reactions with the most variable flux across the chosen strategies."""
        states = [self._strategies.get(n) for n in names]
        states = [s for s in states if s is not None]
        if not states:
            return []
        ids = set()
        for s in states:
            ids |= {k for k, v in s.fluxes.items() if abs(v) > 1e-6}
        def variance(rid):
            vals = [s.flux(rid) for s in states]
            m = sum(vals) / len(vals)
            return sum((v - m) ** 2 for v in vals)
        return sorted(ids, key=variance, reverse=True)[:limit]

    def _render(self) -> None:
        if self._model is None or len(self._strategies) == 0:
            self._placeholder()
            return
        view = self.view_combo.currentIndex()
        if view == 0:
            self._render_difference()
        elif view == 1:
            self._render_single_map()
        elif view == 2:
            self._render_heatmap()
        elif view == 3:
            self._render_waterfall()
        elif view == 4:
            self._render_parallel()

    def _show_flux_map(self, fluxes: dict) -> None:
        ids = [rid for rid, v in fluxes.items() if abs(v) > 1e-9
               and self._model.reactions.has_id(rid)]
        labels = {rid: f"{fluxes[rid]:+.3g}" for rid in ids}
        from ..dialogs.flux_map_dialog import _FluxCategory
        cat = _FluxCategory(name="Flux", reaction_ids=ids)
        self.network.set_categories([cat])
        self.network.set_fluxes(fluxes)
        self.network.set_flux_values(labels)
        self.network.currency_check.setChecked(True)
        self.network.radius_spin.setValue(0)
        self.network.focus_category("Flux")
        self.stack.setCurrentWidget(self.network)

    def _render_difference(self) -> None:
        a, b = self.combo_a.currentText(), self.combo_b.currentText()
        if not a or not b or a == b:
            self.plot.render(plots._empty, "Pick two different strategies (A and B) to see "
                             "the Δflux difference map.", title="difference_map")
            self.stack.setCurrentWidget(self.plot)
            return
        self._show_flux_map(self._strategies.difference(a, b))

    def _render_single_map(self) -> None:
        s = self._strategies.get(self.combo_a.currentText())
        if s is None:
            self._placeholder()
            return
        self._show_flux_map(dict(s.fluxes))

    def _render_heatmap(self) -> None:
        names = self._checked()
        ids = self._top_reactions(names, self.max_rxn.value())
        matrix = self._strategies.flux_matrix(ids)
        self.plot.render(plots.multi_strategy_heatmap, matrix, names, self._label_map(ids),
                         title="multi_strategy_heatmap")
        self.stack.setCurrentWidget(self.plot)

    def _render_waterfall(self) -> None:
        names = self._checked()
        states = [self._strategies.get(n) for n in names]
        rid = self._selected_reaction()
        vals, target = [], ""
        if rid:
            # Follow one chosen reaction across every strategy — this is what makes a
            # heterologous product (absent from the exchange-only target list) plottable.
            for s in states:
                v = s.flux(rid) if s else float("nan")
                vals.append(0.0 if v != v else abs(v))
            target = rid
        else:
            for s in states:
                tf = s.target_flux() if s else float("nan")
                vals.append(0.0 if tf != tf else tf)
                target = target or (s.target if s else "")
        self.plot.render(plots.titre_waterfall, names, vals, target, title="titre_waterfall")
        self.stack.setCurrentWidget(self.plot)

    def _render_parallel(self) -> None:
        names = self._checked()
        ids = self._top_reactions(names, self.max_rxn.value())
        matrix = self._strategies.flux_matrix(ids)
        self.plot.render(plots.parallel_coordinates, matrix, names, self._label_map(ids),
                         title="parallel_coordinates")
        self.stack.setCurrentWidget(self.plot)
