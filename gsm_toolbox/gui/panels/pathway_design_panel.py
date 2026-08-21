"""The 'Pathway Design' tab: predict and apply heterologous pathways.

Workflow: load/fetch one or more reaction databases (managed via a dedicated
dialog, toggled on/off for the search), pick a target metabolite, predict the
heterologous reactions that let the host produce it, then apply that pathway.
The results panel shows one pathway at a time (with a selector when several
alternatives were found) and a row of actions: add it to the model, display it
as a network graph, or explore further alternatives.
"""

from __future__ import annotations

from typing import List

import cobra
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QCompleter,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QPushButton,
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .. import style
from ..views.results_view import ResultsView


class PathwayDesignPanel(QWidget):
    manage_databases_requested = Signal()
    load_database_requested = Signal()
    fetch_database_requested = Signal()
    predict_requested = Signal(str, float)   # target metabolite id, min production
    retrorules_requested = Signal(str)       # target metabolite id — rule-based suggestion
    apply_requested = Signal()
    display_requested = Signal()
    draw_scheme_requested = Signal()
    explore_requested = Signal()
    save_requested = Signal()
    remove_requested = Signal()
    save_design_requested = Signal()
    load_design_requested = Signal()
    database_toggled = Signal(str, bool)   # database name, selected
    load_selected_requested = Signal()
    merge_databases_requested = Signal()
    save_database_requested = Signal()
    balance_steps_requested = Signal()             # balance flagged unbalanced steps
    diagnose_flux_requested = Signal()             # explain why a route carries no flux
    retry_flux_requested = Signal()                # search again for a route THAT RUNS
    fetch_gap_requested = Signal(list)             # fetch the missing chemistry online
    branching_requested = Signal()                 # EA-MNE competition/yield-loss analysis
    mdf_requested = Signal()                       # thermodynamic feasibility (MDF)
    strategies_requested = Signal()                # FSEOF over/under-expression targets
    declare_step_requested = Signal()              # add a final enzyme step to the target
    feasibility_requested = Signal()               # full feasibility report (VI.10)
    structural_search_requested = Signal()         # find targets by structure (VI.2)
    upstream_requested = Signal()                  # explore chemistry feeding the route
    row_context_requested = Signal(list, object)   # forwarded from the active result tab

    def __init__(self):
        super().__init__()
        self._has_db = False
        self._results: List = []
        self._index = 0

        intro = QLabel(
            "Load reaction databases, choose a target metabolite, and predict the "
            "reactions that let the host produce it. Metabolites are matched across "
            "identifier systems; right-click a result for its EC numbers.")
        intro.setWordWrap(True)

        # The three database-management buttons live at the TOP of the right-hand
        # "Reaction databases" panel (added there below) — that keeps them beside the
        # database list they act on, and frees the top row for the target selector.
        self.manage_btn = QPushButton("Manage databases")
        self.manage_btn.setToolTip("See the loaded databases, choose which to include in the "
                                   "search, or permanently delete them.")
        self.manage_btn.clicked.connect(self.manage_databases_requested)
        self.db_btn = QPushButton("Load database")
        self.db_btn.setToolTip("Load a reaction database from a local .json file.")
        self.db_btn.clicked.connect(self.load_database_requested)
        self.fetch_btn = QPushButton("Fetch online")
        self.fetch_btn.setToolTip("Download a reaction database online: BiGG universal or a "
                                  "specific BiGG model, the MetaNetX universal, or a focused "
                                  "KEGG database around a target product.")
        self.fetch_btn.clicked.connect(self.fetch_database_requested)
        from ..widgets.flow_layout import FlowLayout
        db_btn_row = FlowLayout()
        db_btn_row.addWidget(self.manage_btn)
        db_btn_row.addWidget(self.db_btn)
        db_btn_row.addWidget(self.fetch_btn)

        self.target_combo = QComboBox()
        # Wider than before (the database buttons no longer share this row), but bounded:
        # in the wrapping row below, the row's minimum is just the widest single item, so
        # a 280 px combo cannot force the window past the screen.
        self.target_combo.setMinimumWidth(280)
        self.target_combo.setMaximumWidth(460)
        self.target_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # A genome-scale universal adds tens of thousands of metabolites; a
        # contains-match popup completer + select-all-on-focus keeps typing to
        # FILTER the list rather than jump to a prefix or append (#2).
        from ..widgets.dialog_util import configure_search_combo
        configure_search_combo(self.target_combo)
        self.min_flux = QDoubleSpinBox()
        self.min_flux.setRange(0.001, 1000.0)
        self.min_flux.setDecimals(3)
        self.min_flux.setValue(0.1)
        self.min_flux.setToolTip("Minimum target production flux the pathway should enable.")

        self.predict_btn = QPushButton("Predict pathway")
        self.predict_btn.setObjectName("primary")
        self.predict_btn.clicked.connect(self._emit_predict)

        # Rule-based retrosynthesis: propose steps from reaction RULES rather than a
        # database, so it can suggest chemistry no loaded database contains (the answer
        # to targets like prodigiosin). Suggestions are hypotheses to verify.
        self.retrorules_btn = QPushButton("RetroRules Prediction")
        self.retrorules_btn.setToolTip(
            "Suggest heterologous steps from reaction RULES (RetroRules) instead of a "
            "database — it can propose chemistry no database contains. The steps are "
            "predictions to verify, not database-backed reactions.")
        self.retrorules_btn.clicked.connect(self._emit_retrorules)
        # A small second way in, for when the exact name is not in the database (VI.2).
        # The main selector stays literal and predictable; this is the explicit escape.
        self.structural_btn = QPushButton("Structural search…")
        self.structural_btn.setToolTip(
            "Can't find your compound by name? Search the loaded databases by compound "
            "name, SMILES or InChIKey and see the structurally related compounds that "
            "ARE present, with how many reactions can make each.")
        self.structural_btn.clicked.connect(self.structural_search_requested)

        # A wrapping row so the row's minimum is just the widest single control — it can
        # never force the window past the screen, no matter how many buttons it holds.
        target_row = FlowLayout()
        target_row.addWidget(QLabel("Target metabolite:"))
        target_row.addWidget(self.target_combo)
        target_row.addWidget(self.structural_btn)
        target_row.addWidget(QLabel("Min flux:"))
        target_row.addWidget(self.min_flux)
        target_row.addWidget(self.predict_btn)
        target_row.addWidget(self.retrorules_btn)

        self.result_note = QLabel("")
        self.result_note.setWordWrap(True)
        self.result_note.setStyleSheet("font-weight: 600;")
        self.balance_steps_btn = QPushButton("Balance H⁺/H₂O in reactions…")
        self.balance_steps_btn.setToolTip("Add H⁺/H₂O to close protonation/hydration "
                                          "imbalances in the flagged reactions (you confirm each).")
        self.balance_steps_btn.clicked.connect(self.balance_steps_requested)
        self.balance_steps_btn.setVisible(False)
        # Shown only when the selected route cannot carry flux. A zero-flux route is
        # still offered (it is real chemistry, and often the only option) — these say
        # WHY it is zero and offer a way forward, rather than hiding the result.
        self.diagnose_btn = QPushButton("Flux && yield analysis…")
        self.diagnose_btn.setToolTip(
            "Analyse this route: whether a step is blocked (and which one), whether the "
            "objective simply does not reward making the product, and which competing "
            "reactions divert its intermediates away.")
        self.diagnose_btn.clicked.connect(self.diagnose_flux_requested)
        self.diagnose_btn.setVisible(False)
        self.retry_btn = QPushButton("Find a route that runs…")
        self.retry_btn.setToolTip("Search again for an alternative route that can "
                                  "actually carry production flux.")
        self.retry_btn.clicked.connect(self.retry_flux_requested)
        self.retry_btn.setVisible(False)
        # Shown when the search names a compound no loaded database can make. That gap
        # is a DATA problem, so the only real fix is to go and get the chemistry —
        # this fetches it around the missing compound rather than the target.
        self.fetch_gap_btn = QPushButton("Fetch missing chemistry…")
        self.fetch_gap_btn.setToolTip(
            "No loaded database can make one of the precursors this route needs. "
            "Fetch reactions around that compound from KEGG and search again.")
        self.fetch_gap_btn.clicked.connect(self._emit_fetch_gap)
        self.fetch_gap_btn.setVisible(False)
        # Two second-pass checks offered once a concrete route exists. Both are optional:
        # a route is shown and can be added regardless of what they report.
        self.branching_btn = QPushButton("Branching && competition…")
        self.branching_btn.setToolTip(
            "EA-MNE analysis: which native reactions compete for this route's "
            "intermediates, and how much of the ideal yield that diversion costs. "
            "The competitors are candidate knockdown targets.")
        self.branching_btn.clicked.connect(self.branching_requested)
        self.branching_btn.setVisible(False)
        self.mdf_btn = QPushButton("Thermodynamics (MDF)…")
        self.mdf_btn.setToolTip(
            "Estimate whether this route is thermodynamically feasible: per-reaction "
            "ΔrG′ and the pathway Max-min Driving Force (MDF).")
        self.mdf_btn.clicked.connect(self.mdf_requested)
        self.mdf_btn.setVisible(False)
        # Strategy scan works on the HOST-INTEGRATED model, so it is available for
        # rule-based routes too, not only database ones (L9).
        self.strategies_btn = QPushButton("Find strategies (FSEOF)…")
        self.strategies_btn.setToolTip(
            "Scan for engineering targets for this route: which reactions to "
            "overexpress and which to down-regulate to push flux toward the product. "
            "Works for rule-based routes as well as database routes.")
        self.strategies_btn.clicked.connect(self.strategies_requested)
        self.strategies_btn.setVisible(False)
        # Very often a route reaches a compound one known enzyme short of the real
        # target (N-methyltryptamine → DMT). Rather than reporting "not found" for the
        # target, let the user declare that final step and evaluate the whole route (L6).
        self.declare_step_btn = QPushButton("Add a final step…")
        self.declare_step_btn.setToolTip(
            "Your target is often one known enzyme beyond what the database contains. "
            "Declare that final reaction (this route's product → your target) so the "
            "complete pathway can be built and analysed.")
        self.declare_step_btn.clicked.connect(self.declare_step_requested)
        self.declare_step_btn.setVisible(False)
        # The single most important button on a result: everything known about whether
        # this route is actually buildable, in one place (VI.10).
        self.feasibility_btn = QPushButton("Feasibility information…")
        self.feasibility_btn.setToolTip(
            "Full report: chemistry (isomer + balance checks), thermodynamics, flux, "
            "blocking steps and competing reactions.")
        self.feasibility_btn.clicked.connect(self.feasibility_requested)
        self.feasibility_btn.setVisible(False)
        # Alternatives arrive as separate tabs, which makes them hard to weigh against
        # each other. One table, ranked on carbon yield (VI.9).
        self.compare_btn = QPushButton("Compare routes…")
        self.compare_btn.setToolTip(
            "Put every route currently open side by side — steps, carbon yield, flux, "
            "chemistry checks and verdict — ranked by carbon yield.")
        self.compare_btn.clicked.connect(self._compare_routes)
        self.compare_btn.setVisible(False)
        # One core pathway often needs several different terminal steps. Copying it into
        # its own tab lets each variant be edited and analysed without disturbing the rest.
        self.duplicate_btn = QPushButton("Duplicate suggestion")
        self.duplicate_btn.setToolTip(
            "Open an independent copy of this pathway in a new tab, so you can add a "
            "different final step (or edit any step) without changing the original.")
        self.duplicate_btn.clicked.connect(self._duplicate_result)
        self.duplicate_btn.setVisible(False)
        # A route's hand-over to native metabolism is where yield is usually decided, and
        # it is the part the result table says least about — especially when the entry
        # compound turns out to be idle in the host.
        self.upstream_btn = QPushButton("Explore upstream…")
        self.upstream_btn.setToolTip(
            "Look at where this route draws on your host, and find heterologous "
            "reactions that could supply that precursor — essential when the entry "
            "compound carries no flux, and useful for boosting one that does.")
        self.upstream_btn.clicked.connect(self.upstream_requested)
        self.upstream_btn.setVisible(False)
        # The thermodynamics suite is opt-in (Settings ▸ Preferences ▸ Enable MDF Suite).
        # Until then the button never appears, so the feature cannot confuse a user who
        # has not consciously chosen it. The main window sets this on startup and after
        # the preferences dialog closes.
        self._mdf_enabled = False
        # result_note on its own line; the action buttons wrap in a FlowLayout so the
        # row's minimum width is the widest single button, never their sum (the recurring
        # window-width bug: a QHBoxLayout of buttons SUMS its children's minimums).
        note_row = QHBoxLayout()
        note_row.addWidget(self.result_note, 1)
        # The verdict sentence sits directly under the note, before any button, so it is
        # the first thing read about a result (VI.10).
        self.verdict_label = QLabel("")
        self.verdict_label.setWordWrap(True)
        self.verdict_label.setTextFormat(Qt.RichText)
        self.verdict_label.setVisible(False)
        action_btn_row = FlowLayout()
        action_btn_row.addWidget(self.feasibility_btn)
        action_btn_row.addWidget(self.compare_btn)
        action_btn_row.addWidget(self.duplicate_btn)
        action_btn_row.addWidget(self.upstream_btn)
        action_btn_row.addWidget(self.fetch_gap_btn)
        action_btn_row.addWidget(self.diagnose_btn)
        action_btn_row.addWidget(self.branching_btn)
        action_btn_row.addWidget(self.strategies_btn)
        action_btn_row.addWidget(self.declare_step_btn)
        action_btn_row.addWidget(self.mdf_btn)
        action_btn_row.addWidget(self.retry_btn)
        action_btn_row.addWidget(self.balance_steps_btn)

        # One closable tab per analysed/loaded pathway, so several searches stay
        # side by side. Each tab holds a ResultsView of that pathway's reactions.
        self.result_tabs = QTabWidget()
        self.result_tabs.setTabsClosable(True)
        self.result_tabs.setMovable(True)
        self.result_tabs.setDocumentMode(True)
        self.result_tabs.tabCloseRequested.connect(self._close_tab)
        self.result_tabs.currentChanged.connect(lambda _i: self._update_actions_for_current())
        self._tab_views: List[ResultsView] = []

        # bottom action row (short labels with a leading glyph; full text in tooltips)
        self.apply_btn = QPushButton("➕  Add pathway")
        self.apply_btn.setObjectName("primary")
        self.apply_btn.setToolTip("Add this pathway's reactions to the model.")
        self.apply_btn.clicked.connect(self.apply_requested)
        self.display_btn = QPushButton("🕸  Display network")
        self.display_btn.setToolTip("Show this pathway as a network graph, including the native "
                                    "metabolites it connects to.")
        self.display_btn.clicked.connect(self.display_requested)
        self.scheme_btn = QPushButton("🎨  Draw scheme")
        self.scheme_btn.setToolTip("Draw a classic, saveable metabolic-pathway diagram with a "
                                   "structure image for each main metabolite.")
        self.scheme_btn.clicked.connect(self.draw_scheme_requested)
        self.explore_btn = QPushButton("🔍  Alternative pathway")
        self.explore_btn.setToolTip("Omit the pathways shown so far and search for a different "
                                    "route to the same target.")
        self.explore_btn.clicked.connect(self.explore_requested)
        self.remove_btn = QPushButton("🗑  Remove pathway")
        self.remove_btn.setToolTip("Remove the reactions that were added to the model for this "
                                   "pathway.")
        self.remove_btn.clicked.connect(self.remove_requested)
        self.save_btn = QPushButton("💾  Save Results…")
        self.save_btn.setToolTip("Save the predicted pathways to your results folder as a CSV.")
        self.save_btn.clicked.connect(self.save_requested)
        self.save_design_btn = QPushButton("💾  Save")
        self.save_design_btn.setToolTip("Save the designed pathway to a file to reload later.")
        self.save_design_btn.clicked.connect(self.save_design_requested)
        self.load_design_btn = QPushButton("📂  Load")
        self.load_design_btn.setToolTip("Reload a previously saved designed pathway.")
        self.load_design_btn.clicked.connect(self.load_design_requested)
        # Wrapping action rows: never crop, and their minimum width is just the widest
        # single button, so this tab can't force the window wider than the screen.
        from ..widgets.flow_layout import FlowLayout
        action_row = FlowLayout()
        action_row.addWidget(self.apply_btn)
        action_row.addWidget(self.remove_btn)
        action_row.addWidget(self.display_btn)
        action_row.addWidget(self.scheme_btn)
        action_row.addWidget(self.explore_btn)
        action_row2 = FlowLayout()
        action_row2.addWidget(self.save_btn)
        action_row2.addWidget(self.save_design_btn)
        action_row2.addWidget(self.load_design_btn)
        self._set_actions_enabled(False)

        # right-side "available databases" list — tick to select, then Load selected
        self._db_box = QGroupBox("Reaction databases")
        db_v = QVBoxLayout(self._db_box)
        db_v.addLayout(db_btn_row)      # Manage / Load / Fetch, at the top of the panel
        hint = QLabel("Tick the databases you want, then click “Load selected”. "
                      "Ticking does not load anything by itself.")
        hint.setWordWrap(True)
        hint.setStyleSheet(f"color: {style.TEXT_MUTED}; font-size: 11px;")
        db_v.addWidget(hint)
        self._db_list = QListWidget()
        self._db_list.setMaximumWidth(280)
        self._db_list.itemChanged.connect(self._on_db_item_changed)
        db_v.addWidget(self._db_list, 1)
        self.load_selected_btn = QPushButton("Load selected databases")
        self.load_selected_btn.setToolTip("Download (if needed) and load the ticked databases "
                                          "into memory so they can be searched.")
        self.load_selected_btn.clicked.connect(self.load_selected_requested)
        db_v.addWidget(self.load_selected_btn)
        self.merge_btn = QPushButton("Merge loaded databases")
        self.merge_btn.setToolTip("Unify all loaded databases into ONE database where each "
                                  "compound and reaction exists exactly once (duplicates across "
                                  "databases and compartments are removed, collapsed to cytosol).")
        self.merge_btn.clicked.connect(self.merge_databases_requested)
        db_v.addWidget(self.merge_btn)
        self.save_db_btn = QPushButton("Save loaded database…")
        self.save_db_btn.setToolTip("Save the currently loaded (or merged) database to a file so "
                                    "you can reload it offline next time via “Load reaction "
                                    "database…”.")
        self.save_db_btn.clicked.connect(self.save_database_requested)
        db_v.addWidget(self.save_db_btn)
        self._db_totals = QLabel("No databases loaded.")
        self._db_totals.setWordWrap(True)
        self._db_totals.setStyleSheet(f"color: {style.TEXT_MUTED};")
        db_v.addWidget(self._db_totals)

        center = QHBoxLayout()
        left = QVBoxLayout()
        left.addLayout(note_row)
        left.addWidget(self.verdict_label)
        left.addLayout(action_btn_row)
        left.addWidget(self.result_tabs, 1)
        left.addLayout(action_row)
        left.addLayout(action_row2)
        center.addLayout(left, 1)
        center.addWidget(self._db_box)

        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(target_row)
        layout.addLayout(center, 1)
        self._set_enabled(False)

    # ----- enable state ------------------------------------------------
    def _set_enabled(self, on: bool) -> None:
        self.target_combo.setEnabled(on)
        self.min_flux.setEnabled(on)
        self.predict_btn.setEnabled(on)
        self.retrorules_btn.setEnabled(on)

    def _set_actions_enabled(self, on: bool) -> None:
        self.apply_btn.setEnabled(on)
        self.display_btn.setEnabled(on)
        self.scheme_btn.setEnabled(on)
        self.explore_btn.setEnabled(on)
        self.save_btn.setEnabled(on)
        self.remove_btn.setEnabled(on)
        self.save_design_btn.setEnabled(on)
        # "Load designed pathway" is always available (even with nothing shown yet).
        self.load_design_btn.setEnabled(True)

    def set_model(self, model: cobra.Model) -> None:
        # Target choices are supplied by the main window via set_targets(). Do NOT
        # wipe the predicted-pathway tabs here: this runs on every model refresh
        # (including right after "Add pathway to model"), and the user must keep
        # interacting with the displayed pathways (display/scheme/explore/remove).
        self._update_actions_for_current()

    def clear_results(self) -> None:
        """Discard all pathway tabs (used when a brand-new model is opened)."""
        self.result_tabs.clear()
        self._results = []
        self._tab_views = []
        self.result_note.setText("")
        self._set_actions_enabled(False)

    def set_databases_summary(self, entries: list) -> None:
        """Show the available databases as checkboxes.

        ``entries`` is ``[{name, reactions, metabolites, selected, loaded, source}]``.
        Ticking a row selects it (emits ``database_toggled``) but does not load it —
        loading happens when the user clicks “Load selected databases”. A loaded DB
        shows its real reaction count; an unloaded one shows an approximate size."""
        from PySide6.QtWidgets import QListWidgetItem
        self._db_list.blockSignals(True)
        self._db_list.clear()
        n_loaded = n_selected = 0
        for e in entries:
            if e.get("loaded"):
                detail = f"{e['reactions']:,} rxns · loaded"
            else:
                detail = f"~{e['reactions']:,} rxns · not loaded"
            item = QListWidgetItem(f"{e['name']}   ({detail})")
            item.setFlags((item.flags() | Qt.ItemIsUserCheckable) & ~Qt.ItemIsSelectable)
            item.setCheckState(Qt.Checked if e.get("selected") else Qt.Unchecked)
            item.setData(Qt.UserRole, e["name"])
            self._db_list.addItem(item)
            n_selected += bool(e.get("selected"))
            n_loaded += bool(e.get("loaded"))
        self._db_list.blockSignals(False)
        # The search can run once at least one selected DB is loaded.
        self._has_db = any(e.get("selected") and e.get("loaded") for e in entries)
        if not entries:
            self._db_totals.setText("No databases available.")
        else:
            self._db_totals.setText(
                f"{n_selected} selected · {n_loaded} loaded of {len(entries)} available.")
        self._set_enabled(self._has_db)

    def _on_db_item_changed(self, item) -> None:
        name = item.data(Qt.UserRole)
        if name:
            self.database_toggled.emit(str(name), item.checkState() == Qt.Checked)

    def set_targets(self, items: list) -> None:
        """Populate the target list. ``items`` is a list of (display, id) pairs."""
        current = self.target_combo.currentData()
        self._target_ids = {item_id for _, item_id in items}
        self.target_combo.blockSignals(True)
        self.target_combo.clear()
        for display, item_id in items:
            self.target_combo.addItem(display, item_id)
        if current is not None:
            idx = self.target_combo.findData(current)
            if idx >= 0:
                self.target_combo.setCurrentIndex(idx)
        self.target_combo.blockSignals(False)

    def _resolve_target(self) -> str:
        idx = self.target_combo.currentIndex()
        text = self.target_combo.currentText().strip()
        if idx >= 0 and self.target_combo.itemText(idx) == text:
            data = self.target_combo.itemData(idx)
            if data:
                return str(data)
        if text in getattr(self, "_target_ids", set()):
            return text
        if text.endswith(")") and "(" in text:
            return text[text.rfind("(") + 1:-1].strip()
        return text

    def current_target(self) -> str:
        return self._resolve_target()

    def _emit_predict(self) -> None:
        target = self._resolve_target()
        if target:
            self.predict_requested.emit(target, self.min_flux.value())

    def _emit_retrorules(self) -> None:
        target = self._resolve_target()
        if target:
            self.retrorules_requested.emit(target)

    # ----- results (one closable tab per pathway) ----------------------
    def add_results(self, results: List) -> None:
        """Append a tab for each non-empty pathway (older tabs are kept)."""
        added = 0
        for r in results or []:
            if r is None or r.reactions.empty:
                continue
            self._add_tab(r)
            added += 1
        if added == 0 and not self._results:
            # Nothing found and nothing shown: reflect the note (no tab to add).
            nf = next((r for r in (results or []) if r is not None), None)
            self.result_note.setText(nf.note if nf else "No pathway found.")
            # A not-found result has no tab, so the tab-driven action update never
            # runs for it — yet this is exactly when the fetch-the-missing-chemistry
            # offer matters most. Drive it from the result itself.
            self._last_not_found = nf
            self.fetch_gap_btn.setVisible(bool(getattr(nf, "missing_compounds", None)))
            self.diagnose_btn.setVisible(False)
            self.retry_btn.setVisible(False)
            self.balance_steps_btn.setVisible(False)
            return
        if added:
            self._last_not_found = None
            self.result_tabs.setCurrentIndex(self.result_tabs.count() - 1)
        self._update_actions_for_current()

    def append_result(self, result) -> None:
        """Add a newly-explored alternative as its own tab and select it."""
        if result is not None and not result.reactions.empty:
            self._add_tab(result)
            self.result_tabs.setCurrentIndex(self.result_tabs.count() - 1)
        self._update_actions_for_current()

    def _tab_title(self, result) -> str:
        base = result.target.rsplit("_", 1)[0] if "_" in result.target else result.target
        n = len(result.reaction_ids)
        # Disambiguate multiple pathways to the same target (#, #2, #3…).
        existing = sum(1 for r in self._results
                       if r.target == result.target)
        suffix = f" #{existing + 1}" if existing else ""
        return f"{base}{suffix} ({n})"

    def _add_tab(self, result, *, copied: bool = False) -> None:
        view = ResultsView()
        view.row_context_requested.connect(self.row_context_requested)
        view.show_dataframe(result.reactions,
                            f"Heterologous reactions to add for {result.target}")
        # The suggested id is user-editable — the name the reaction gets on Add (#B5).
        view.set_editable_columns({"suggested_id"})
        title = self._tab_title(result)     # count same-target tabs BEFORE appending
        if copied:
            # A copy is a variant of an existing route, not a newly found alternative;
            # marking it keeps the two apart in the tab bar and in Compare routes.
            title += " ✎"
        self._results.append(result)
        self._tab_views.append(view)
        self.result_tabs.addTab(view, title)

    def _close_tab(self, index: int) -> None:
        if not (0 <= index < len(self._results)):
            return
        self.result_tabs.removeTab(index)
        del self._results[index]
        del self._tab_views[index]
        self._update_actions_for_current()

    def _update_actions_for_current(self) -> None:
        result = self.current_result()
        if result is None:
            self.result_note.setText("")
            self.compare_btn.setVisible(False)
            self.duplicate_btn.setVisible(False)
            self.upstream_btn.setVisible(False)
            self._set_actions_enabled(False)
            return
        f = result.production_flux
        has_flux = f == f          # not NaN
        if has_flux and getattr(result, "flux_is_indicative", False):
            # Rule-derived route: report only whether it can carry flux. The absolute
            # number is an artefact of a synthetic demand on a novel compound (L2).
            flux_note = ("  This rule-based route <b>can carry flux</b> (the absolute "
                         "value is not a validated capacity — compare rule routes only "
                         "to each other)." if f > 1e-9 else
                         "  This rule-based route <b>carries no flux</b> as written.")
        elif has_flux:
            flux_note = f"  Predicted production flux: {f:.4g}."
            cy = getattr(result, "carbon_yield", float("nan"))
            if cy == cy:           # normalised, cross-target-comparable figure (L10)
                flux_note += f"  Carbon yield: {cy * 100:.1f}% of consumed C."
            elif getattr(result, "carbon_yield_note", ""):
                # Say why there is no yield. A blank (or a NaN) reads as a failed design;
                # "the product has no formula" is a database gap the user can close.
                flux_note += ("  Carbon yield: not computable — "
                              f"{result.carbon_yield_note}.")
        else:
            flux_note = ""
        self.result_note.setText(f"Target {result.target}: {result.note}{flux_note}")
        self.result_note.setTextFormat(Qt.RichText)
        self._set_actions_enabled(True)
        # Offer the balance shortcut only when the route has unbalanced/unverified steps.
        self.balance_steps_btn.setVisible(not getattr(result, "balanced", True))
        # A route with steps but no predicted flux is still offered — explain it and
        # offer a way on, instead of leaving the user with an unexplained 0.
        f = result.production_flux
        no_flux = bool(result.reaction_ids) and not (f == f and f > 1e-9)
        # Flux/yield analysis is useful for ANY route (a working route still has
        # branches competing for its intermediates); the retry only makes sense when
        # the route cannot run.
        self.diagnose_btn.setVisible(bool(result.reaction_ids))
        # Branching is a second-pass check on a concrete route: offer it whenever there
        # are heterologous steps. Thermodynamics is additionally gated on the opt-in
        # preference, so it stays invisible unless the user enabled the MDF suite.
        self.branching_btn.setVisible(bool(result.reaction_ids))
        self.strategies_btn.setVisible(bool(result.reaction_ids))
        self.declare_step_btn.setVisible(bool(result.reaction_ids))
        self.feasibility_btn.setVisible(bool(result.reaction_ids))
        # Comparison only says anything with at least two routes on the table.
        self.compare_btn.setVisible(len(self._results) > 1)
        self.duplicate_btn.setVisible(bool(result.reaction_ids))
        self.upstream_btn.setVisible(bool(result.reaction_ids))
        self._update_verdict(result)
        self.mdf_btn.setVisible(bool(result.reaction_ids) and self._mdf_enabled)
        self.retry_btn.setVisible(no_flux)
        self.fetch_gap_btn.setVisible(bool(getattr(result, "missing_compounds", None)))

    def _duplicate_result(self) -> None:
        """Copy the selected route into its own tab so it can be edited independently."""
        result = self.current_result()
        if result is None or not result.reaction_ids:
            return
        clone = result.duplicate()
        clone.note = f"Copy of “{self._tab_title(result)}”. {clone.note}"
        self._add_tab(clone, copied=True)
        self.result_tabs.setCurrentIndex(self.result_tabs.count() - 1)
        self._update_actions_for_current()

    def _compare_routes(self) -> None:
        """Open the side-by-side route table (VI.9). Read-only, so it needs no model."""
        if len(self._results) < 2:
            return
        from ..dialogs.compare_routes_dialog import CompareRoutesDialog
        titles = [self.result_tabs.tabText(i) for i in range(len(self._results))]
        CompareRoutesDialog(self, list(self._results), titles).exec()

    def _update_verdict(self, result) -> None:
        """Show the one-sentence buildability verdict for the selected route (VI.10).

        Only what the search itself computed is used here (chemistry, balance, flux) —
        thermodynamics and competition are added once the user runs those analyses, and
        the sentence says plainly when they are not yet known.
        """
        if result is None or not getattr(result, "reaction_ids", None):
            self.verdict_label.setVisible(False)
            return
        try:
            from ...core import feasibility as fz
            rep = fz.assess(result,
                            diagnosis=getattr(result, "_diagnosis", None),
                            branching=getattr(result, "_branching", None),
                            mdf=getattr(result, "_mdf", None))
            self.verdict_label.setText(
                f"<span style='color:{rep.colour}'><b>{rep.label}</b></span> — "
                + rep.sentence().split("— ", 1)[-1])
            self.verdict_label.setVisible(True)
        except Exception:  # noqa: BLE001 — a verdict must never break the panel
            self.verdict_label.setVisible(False)

    def set_mdf_enabled(self, on: bool) -> None:
        """Show/hide the thermodynamics controls (Settings ▸ Preferences opt-in)."""
        self._mdf_enabled = bool(on)
        self._update_actions_for_current()

    def refresh_action_visibility(self) -> None:
        """Re-evaluate which result-action buttons are shown (after a preference change)."""
        self._update_actions_for_current()

    def _emit_fetch_gap(self) -> None:
        # The gap normally belongs to a not-found result, which has no tab of its own.
        result = self.current_result() or getattr(self, "_last_not_found", None)
        missing = list(getattr(result, "missing_compounds", None) or []) if result else []
        if missing:
            self.fetch_gap_requested.emit(missing)

    def current_result(self):
        i = self.result_tabs.currentIndex()
        if 0 <= i < len(self._results):
            return self._results[i]
        return None

    def refresh_current_result(self) -> None:
        """Re-render the current tab's table from its (possibly mutated) result —
        used after balancing a suggested reaction in place."""
        i = self.result_tabs.currentIndex()
        if 0 <= i < len(self._tab_views):
            r = self._results[i]
            self._tab_views[i].show_dataframe(
                r.reactions, f"Heterologous reactions to add for {r.target}")
            self._tab_views[i].set_editable_columns({"suggested_id"})
        self._update_actions_for_current()

    def all_reaction_ids(self) -> set:
        ids: set = set()
        for r in self._results:
            ids |= set(r.reaction_ids)
        return ids

    def all_results(self) -> list:
        return list(self._results)

    # Entry points used by the main window.
    def show_results(self, results: List) -> None:
        self.add_results(results)

    def show_result(self, result) -> None:
        self.add_results([result])


def _empty_df():
    import pandas as pd
    return pd.DataFrame()
