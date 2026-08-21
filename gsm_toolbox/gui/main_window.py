"""The main application window: menus, docks, central tabs, and analysis wiring.

Layout:
* Center  — tabbed workspace: Explorer (reaction/metabolite/gene lists, the main
  panel), Network Map, Analysis.
* Left    — Categories panel.
* Right   — Info panel: a Model tab (general info) and a Selection tab (details of
  the clicked reaction/metabolite/gene).

Analyses live only in the Analysis tab (not the menu bar) and always run on a
private copy of the model on a worker thread, which keeps the UI responsive and
avoids GLPK's non-thread-safe shared-solver crashes.
"""

from __future__ import annotations

import os
from contextlib import contextmanager
from typing import Callable, Optional

import cobra
from PySide6.QtCore import QEvent, QSize, Qt, QTimer
from PySide6.QtGui import QAction, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMenu,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .. import __app_name__, __version__
from ..core import categories as cat_core
from ..core import community as community_core
from ..core import (cache, databases, editing, io_models, media, namespace,
                    pathway_design, pathway_search)
from ..core.analysis import deletions, fba, gapfill, mutant, omics, pathways, phenotype, qc
from ..core.analysis import strain_design as sd
from ..core.project import PROJECT_EXT, Project, ProjectError
from .dialogs.media_dialog import MediaDialog
from .dialogs.reaction_builder import ReactionBuilderDialog
from .dialogs.reaction_dialog import ReactionDialog
from .panels.analysis_panel import AnalysisPanel
from .panels.categories_panel import CategoriesPanel
from .panels.growth_settings_panel import GrowthSettingsPanel
from .panels.info_panel import InfoPanel
from .panels.model_explorer import ModelExplorer
from .panels.pathway_design_panel import PathwayDesignPanel
from .panels.settings_panel import SettingsPanel
from .panels.strategy_explorer import StrategyExplorer
from .panels.escher_explorer import EscherExplorer
from .views.network_view import NetworkView
from .views.results_view import ResultsView
from .widgets.busy import run_busy, was_cancelled
from .widgets.job_status import JobStatusBar
from .workers import JobManager

#: Qt's QWIDGETSIZE_MAX — "no maximum".
_QWIDGETSIZE_MAX = 16777215

_MODEL_FILTER = "Models (*.xml *.sbml *.json *.mat);;All files (*)"
_PROJECT_FILTER = f"GSM ToolBox project (*{PROJECT_EXT})"


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.project: Optional[Project] = None
        self._last_result: Optional[fba.FBAResult] = None
        self.jobs = JobManager(self)
        self._results_dir: Optional[str] = None      # auto-save analysis result tables here
        self.community: Optional[community_core.CommunityModel] = None

        # Window-state guards — see hold_window_state and minimumSizeHint.
        self._layout_ops = 0
        self._settling = False
        self._should_be_maximized = False
        self._user_touched_frame = False

        self.setWindowTitle(__app_name__)
        self.resize(1360, 860)
        # Allow docks to be nested/split freely so floated panels can be dragged
        # back into any position in the main window.
        self.setDockNestingEnabled(True)

        self._build_docks()
        self._build_central()
        self._wire_signals()
        self._build_actions()
        self._build_statusbar()
        self._arrange_docks()
        self._apply_appearance()          # font size, info-panel side, dock visibility
        self._update_action_availability()   # calls _update_action_states in turn

    def resizeEvent(self, event):  # noqa: N802 - Qt override (A1)
        """Keep the whole window inside the screen at its normal size, without ever
        fighting the maximized/fullscreen state (see _cap_to_screen)."""
        super().resizeEvent(event)
        self._cap_to_screen()

    def moveEvent(self, event):  # noqa: N802 - Qt override
        # Moving between monitors changes which screen (and cap) applies.
        super().moveEvent(event)
        self._cap_to_screen()

    #: Events that mean the *user* is acting on the window frame — the restore button,
    #: a title-bar double-click, a drag off the top edge. Qt does not say who caused a
    #: WindowStateChange, so this is how the two are told apart.
    _FRAME_INPUT = {
        QEvent.NonClientAreaMouseButtonPress,
        QEvent.NonClientAreaMouseButtonDblClick,
        QEvent.NonClientAreaMouseButtonRelease,
    }

    def event(self, event):  # noqa: N802 - Qt override
        """Record genuine user interaction with the window frame.

        An earlier attempt guessed from timing — treat a state change during a relayout as
        spurious, anything at idle as the user. That failed because relayouts settle
        several event-loop turns later (a dock animating, a modal progress dialog closing,
        a panel repopulating), by which point the guard had expired and a programmatic
        un-maximize was mistaken for a deliberate one. Watching for input is exact:
        un-maximizing needs a click on the frame or a keystroke, and nothing else can
        produce one.
        """
        if event.type() in self._FRAME_INPUT:
            self._note_frame_input()
        elif event.type() == QEvent.KeyPress:
            # Win+Down and Alt+Space are the keyboard routes out of Maximized. Plain
            # keystrokes are not — typing in a search box must not be read as consent to
            # un-maximize.
            if event.modifiers() & (Qt.MetaModifier | Qt.AltModifier):
                self._note_frame_input()
        return super().event(event)

    def _note_frame_input(self) -> None:
        """Arm the "the user did this" flag, briefly.

        It expires because a click on the frame authorises the state change that follows
        it immediately, not one that happens a minute later during a model load.
        """
        self._user_touched_frame = True
        QTimer.singleShot(400, self._clear_frame_input)

    def _clear_frame_input(self) -> None:
        self._user_touched_frame = False

    def changeEvent(self, event):  # noqa: N802 - Qt override
        super().changeEvent(event)
        if event.type() != QEvent.WindowStateChange:
            return
        if self.isMaximized() or self.isFullScreen():
            # However it got there, maximized is now the state to keep.
            self._should_be_maximized = True
            self._user_touched_frame = False
        elif self._user_touched_frame:
            # The user clicked restore, double-clicked the title bar or used the
            # keyboard. Honour it and stop re-maximizing.
            self._should_be_maximized = False
            self._user_touched_frame = False
        else:
            # Nothing the user did asked for this — a relayout knocked it out. Put it
            # back once the current batch of geometry changes has drained.
            QTimer.singleShot(0, self._restore_maximized_if_lost)
        QTimer.singleShot(0, self._cap_to_screen)

    # ----- keeping the window maximized ---------------------------------
    def minimumSizeHint(self) -> QSize:  # noqa: N802 - Qt override
        """Never report a minimum larger than the screen.

        This is the root cause of the "window un-maximizes by itself" bug. Qt honours a
        layout minimum over the maximized state, so the moment any child's minimum pushes
        the total past the work area — a combo box filled with long reaction names when a
        model loads, a dock reappearing — Qt silently restores the window to fit. Clamping
        here means that situation cannot arise; content scrolls or elides instead.
        """
        hint = super().minimumSizeHint()
        screen = self.screen() or QGuiApplication.primaryScreen()
        if screen is None:
            return hint
        avail = screen.availableGeometry()
        return QSize(min(hint.width(), avail.width()),
                     min(hint.height(), avail.height()))

    @contextmanager
    def hold_window_state(self):
        """Run a relayout without letting it change the maximized state.

        Loading a model, toggling a dock and repopulating a panel all resize widgets, and
        any of them can knock the window out of Maximized. Only the user should be able to
        do that, so programmatic work is bracketed and the state restored afterwards.
        """
        was_maximized = self.isMaximized() or self.isFullScreen()
        self._layout_ops += 1
        try:
            yield
        finally:
            self._layout_ops -= 1
            if was_maximized:
                self._should_be_maximized = True
                self._restore_maximized_if_lost()
                # Some relayouts settle a tick later (dock animations, tab switches), so
                # check again once the event loop has drained.
                self._settling = True
                QTimer.singleShot(0, self._finish_settling)

    def _finish_settling(self) -> None:
        self._restore_maximized_if_lost()
        self._settling = False

    def _restore_maximized_if_lost(self) -> None:
        if self._should_be_maximized and not (self.isMaximized() or self.isFullScreen()):
            # Lift any cap first: a maximized window is a few px larger than the work
            # area on Windows, so a stale maximumSize would bounce it straight back out.
            self.setMaximumSize(_QWIDGETSIZE_MAX, _QWIDGETSIZE_MAX)
            self.showMaximized()

    def _on_dock_visibility_changed(self, *_args) -> None:
        """Closing or reopening a dock must not change the window state."""
        if self._should_be_maximized:
            QTimer.singleShot(0, self._restore_maximized_if_lost)

    def _cap_to_screen(self) -> None:
        """Guarantee the window is fully on-screen at its normal size, WITHOUT ever
        constraining a maximized/fullscreen window (a Windows maximized window is a
        few px larger than the work area, so a maximumSize cap makes it silently
        restore — the "running an analysis / switching tabs drops out of Maximized"
        bug). While maximized we lift the cap; in the normal state we cap the size to
        the work area and nudge the window back on-screen if it spilled past an edge.
        Enforcement is done purely via setMaximumSize/move — never resize() inside the
        resize event — to avoid re-entrancy."""
        _MAXQ = _QWIDGETSIZE_MAX
        if self.isMaximized() or self.isFullScreen():
            if self.maximumWidth() != _MAXQ or self.maximumHeight() != _MAXQ:
                self.setMaximumSize(_MAXQ, _MAXQ)
            return
        if getattr(self, "_clamping", False):
            return
        if self._should_be_maximized:
            # The window is meant to be maximized and something knocked it out. Put it
            # back rather than capping it here — an earlier version returned at this
            # point, which is why the stray restore was left visibly larger than the
            # screen instead of being corrected.
            self._restore_maximized_if_lost()
            return
        screen = self.screen() or QGuiApplication.primaryScreen()
        if screen is None:
            return
        avail = screen.availableGeometry()
        self._clamping = True
        try:
            if self.maximumWidth() != avail.width() or self.maximumHeight() != avail.height():
                self.setMaximumSize(avail.width(), avail.height())
            # Shrink first if the window is genuinely wider or taller than the screen.
            # The cap above stops it *growing*, but a geometry Qt restored from a stale
            # oversized "normal" size is already too big and has to be brought back — the
            # "restored view spanning more than the screen" symptom.
            fg = self.frameGeometry()
            if fg.width() > avail.width() or fg.height() > avail.height():
                frame = fg.size() - self.size()      # title bar and borders
                self.resize(max(avail.width() - frame.width(), self.minimumWidth()),
                            max(avail.height() - frame.height(), self.minimumHeight()))
            fg = self.frameGeometry()
            if fg.width() <= avail.width() and fg.height() <= avail.height() \
                    and not avail.contains(fg):
                x = min(max(fg.x(), avail.x()), avail.right() - fg.width() + 1)
                y = min(max(fg.y(), avail.y()), avail.bottom() - fg.height() + 1)
                if (x, y) != (fg.x(), fg.y()):
                    self.move(max(x, avail.x()), max(y, avail.y()))
        finally:
            self._clamping = False

    # ----- construction ------------------------------------------------
    def _build_central(self) -> None:
        """Center holds the Network Map and Analysis tabs (the main work area)."""
        self.tabs = QTabWidget()
        self.network_view = NetworkView()
        self.analysis_panel = AnalysisPanel()
        self.growth_panel = GrowthSettingsPanel()
        self.pathway_panel = PathwayDesignPanel()
        self.strategy_explorer = StrategyExplorer()
        self.escher_explorer = EscherExplorer()
        from .panels.dynamic_panel import DynamicAnalysisPanel
        self.dynamic_panel = DynamicAnalysisPanel()
        # Four tabs, each a distinct kind of work. Growth Settings moved to the toolbar
        # and the two visualizers behind one Network Visualization button, which is what
        # made room for Dynamic Analysis without the strip becoming unreadable.
        self.tabs.addTab(self.analysis_panel, "Analysis")
        self.tabs.addTab(self.pathway_panel, "Pathway Design")
        self.tabs.addTab(self.dynamic_panel, "Dynamic Analysis")
        self.tabs.addTab(self.network_view, "Network Map")
        self.dynamic_panel.scan_requested.connect(self._run_condition_scan)
        self.dynamic_panel.timecourse_requested.connect(self._run_timecourse)
        self.dynamic_panel.regulation_requested.connect(self._open_regulation)
        self.dynamic_panel.plot_requested.connect(self._plot_dynamic_results)
        self.dynamic_panel.table_requested.connect(self._show_dynamic_run)
        # Clicking a stored run opens it; that is what the tab is for.
        self.dynamic_panel.run_tabs.tabBarClicked.connect(
            lambda i: self._show_dynamic_run(self.dynamic_panel.run_tabs.tabText(i)))
        self.growth_panel.apply_requested.connect(self._apply_growth_settings)
        self.growth_panel.mode_requested.connect(self._apply_growth_mode)
        self.strategy_explorer.save_strategy_requested.connect(self._save_strategy)
        # Run FBA/pFBA from the Strategy tab itself (no tab round-trip just to solve).
        self.strategy_explorer.run_analysis_requested.connect(self._run_analysis)
        self.strategy_explorer.remove_strategy_requested.connect(self._remove_strategy)
        self.escher_explorer.reaction_info_requested.connect(self._show_reaction_details)
        self.escher_explorer.metabolite_info_requested.connect(self._show_metabolite_details)
        self.tabs.setCurrentWidget(self.analysis_panel)   # Analysis is the default view
        self.tabs.currentChanged.connect(self._on_tab_changed)
        self.setCentralWidget(self.tabs)
        self._databases: list = []          # [{name, model, included, cache_path, source}]
        self._pathway_db: Optional[object] = None   # combined model of included databases
        self._pathway_results = []
        self._applied_by_target: dict = {}   # target -> reaction ids added (removal fallback)
        self._flux_data = None                # {reaction_id: (flux, label)} from last flux table
        self._startup_defaults_done = False
        self.omics_panel = None               # the Omics tab, created lazily (#F4)

    def _on_tab_changed(self, index: int) -> None:
        # The first time the Pathway Design tab is opened, just LIST the available
        # reaction databases (unchecked, no loading — instant). The user ticks the
        # ones they want and clicks "Load selected databases"; heavy parsing/merging
        # happens then, and only at prediction time is the combined universal built.
        if self.tabs.widget(index) is not self.pathway_panel:
            return
        if self._startup_defaults_done or os.environ.get("GSM_SELFTEST"):
            return
        self._startup_defaults_done = True
        self._populate_available_databases()

    def _build_docks(self) -> None:
        """Left column = Explorer (top) + Categories (bottom); Information = bottom (compact)."""
        features = (QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable
                    | QDockWidget.DockWidgetClosable)

        self.explorer = ModelExplorer()
        self.explorer_dock = QDockWidget("Explorer", self)
        self.explorer_dock.setObjectName("explorer_dock")
        self.explorer_dock.setWidget(self.explorer)
        self.explorer_dock.setAllowedAreas(Qt.AllDockWidgetAreas)
        self.explorer_dock.setFeatures(features)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.explorer_dock)

        self.categories_panel = CategoriesPanel()
        self.categories_dock = QDockWidget("Categories", self)
        self.categories_dock.setObjectName("categories_dock")
        self.categories_dock.setWidget(self.categories_panel)
        self.categories_dock.setAllowedAreas(Qt.AllDockWidgetAreas)
        self.categories_dock.setFeatures(features)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.categories_dock)
        # Stack categories under the explorer in the left column.
        self.splitDockWidget(self.explorer_dock, self.categories_dock, Qt.Vertical)

        # Information panel: docked on the RIGHT by default (most comfortable for
        # reading details next to the map/tables). It can be collapsed/closed and
        # re-opened from the View menu, and moved anywhere.
        self.info = InfoPanel()
        self.info_dock = QDockWidget("Information", self)
        self.info_dock.setObjectName("info_dock")
        self.info_dock.setWidget(self.info)
        self.info_dock.setAllowedAreas(Qt.AllDockWidgetAreas)
        self.info_dock.setFeatures(features)
        self.addDockWidget(Qt.RightDockWidgetArea, self.info_dock)
        # Hidden by default so the central tables/maps get the full width and the
        # window fits the screen (a Model Summary pops up on load instead). It stays
        # available from View ▸ Information, and clicking anything still fills it.
        self.info_dock.hide()

        # Showing or closing a dock relayouts the whole window, which used to knock it
        # out of Maximized. Watch each one and put the state back.
        for dock in (self.explorer_dock, self.categories_dock, self.info_dock):
            dock.visibilityChanged.connect(self._on_dock_visibility_changed)
            dock.topLevelChanged.connect(self._on_dock_visibility_changed)

    def _apply_appearance(self) -> None:
        """Apply the Appearance preferences to the running window.

        Called once at startup and again whenever Preferences is accepted, so a change
        takes effect without a restart. Every value is read defensively: a preferences
        file edited by hand should never stop the window from opening.
        """
        from ..core import preferences as prefs

        try:
            size = int(prefs.get(prefs.FONT_SIZE) or 9)
        except (TypeError, ValueError):
            size = 9
        size = max(7, min(size, 18))
        font = self.font()
        if font.pointSize() != size:
            font.setPointSize(size)
            self.setFont(font)
            app = QApplication.instance()
            if app is not None:
                app_font = app.font()
                app_font.setPointSize(size)
                app.setFont(app_font)

        where = prefs.get(prefs.INFO_PANEL_POSITION) or "right"
        area = Qt.RightDockWidgetArea if where == "right" else Qt.BottomDockWidgetArea
        if self.dockWidgetArea(self.info_dock) != area:
            self.addDockWidget(area, self.info_dock)
        for action, matches in ((self.act_info_right, where == "right"),
                                (self.act_info_bottom, where != "right")):
            action.setChecked(matches)

        # Dock visibility. The Information dock starts hidden by design, so it is only
        # forced open when the user has explicitly asked for it.
        for dock, key in ((self.explorer_dock, prefs.SHOW_EXPLORER),
                          (self.categories_dock, prefs.SHOW_CATEGORIES),
                          (self.info_dock, prefs.SHOW_INFO)):
            dock.setVisible(bool(prefs.get(key)))

    def _arrange_docks(self) -> None:
        """Give the docks sensible initial sizes: a compact bottom info strip and a
        moderate left column, leaving the center maximally wide for tables/maps."""
        self.resizeDocks([self.explorer_dock, self.categories_dock], [520, 240], Qt.Vertical)
        self.resizeDocks([self.explorer_dock], [320], Qt.Horizontal)
        self.resizeDocks([self.info_dock], [340], Qt.Horizontal)  # comfortable right column

    def _wire_signals(self) -> None:
        # Explorer
        self.explorer.reaction_selected.connect(self._on_reaction_selected)
        self.explorer.metabolite_selected.connect(self.info.show_metabolite)
        self.explorer.metabolite_details_requested.connect(self._show_metabolite_details)
        self.explorer.edit_metabolite_requested.connect(self._edit_metabolite)
        self.explorer.gene_selected.connect(self.info.show_gene)
        self.explorer.edit_reaction_requested.connect(self._edit_reaction)
        self.explorer.reaction_details_requested.connect(self._show_reaction_details)
        self.explorer.show_in_map_requested.connect(self._show_in_map)
        self.explorer.add_to_category_requested.connect(self._add_reactions_to_category)
        # Map
        self.network_view.node_selected.connect(self._on_network_node)
        self.network_view.node_context_requested.connect(self._on_node_context)
        # Analysis
        self.analysis_panel.run_requested.connect(self._run_analysis)
        self.analysis_panel.save_requested.connect(self._save_analysis_results)
        self.analysis_panel.display_fluxes_requested.connect(self._display_fluxes)
        self.analysis_panel.visualize_requested.connect(self._open_plot_gallery)
        self.analysis_panel.objective_bar.objective_changed.connect(self._apply_objective)
        self.analysis_panel.objective_bar.consortia_requested.connect(self._edit_consortia_objective)
        self.analysis_panel.results_view.row_context_requested.connect(self._on_result_context)
        # Pathway design
        self.pathway_panel.load_database_requested.connect(self._load_pathway_database)
        self.pathway_panel.fetch_database_requested.connect(self._fetch_pathway_database)
        self.pathway_panel.predict_requested.connect(self._predict_pathway)
        self.pathway_panel.retrorules_requested.connect(self._retrorules_suggest)
        self.pathway_panel.apply_requested.connect(self._apply_pathway)
        self.pathway_panel.display_requested.connect(self._display_pathway)
        self.pathway_panel.draw_scheme_requested.connect(self._draw_pathway_scheme)
        self.pathway_panel.explore_requested.connect(self._explore_alternative_pathway)
        self.pathway_panel.save_requested.connect(self._save_pathway_results)
        self.pathway_panel.remove_requested.connect(self._remove_pathway)
        self.pathway_panel.save_design_requested.connect(self._save_designed_pathway)
        self.pathway_panel.load_design_requested.connect(self._load_designed_pathway)
        self.pathway_panel.database_toggled.connect(self._set_database_included)
        self.pathway_panel.load_selected_requested.connect(self._load_selected_databases)
        self.pathway_panel.merge_databases_requested.connect(self._merge_databases)
        self.pathway_panel.save_database_requested.connect(self._save_loaded_database)
        self.pathway_panel.balance_steps_requested.connect(self._balance_current_pathway_steps)
        self.pathway_panel.diagnose_flux_requested.connect(self._diagnose_pathway_flux)
        self.pathway_panel.retry_flux_requested.connect(self._retry_for_flux_carrying_route)
        self.pathway_panel.fetch_gap_requested.connect(self._fetch_missing_chemistry)
        self.pathway_panel.branching_requested.connect(self._branching_analysis)
        self.pathway_panel.mdf_requested.connect(self._mdf_analysis)
        self.pathway_panel.strategies_requested.connect(self._pathway_strategies)
        self.pathway_panel.declare_step_requested.connect(self._declare_final_step)
        self.pathway_panel.feasibility_requested.connect(self._feasibility_report)
        self.pathway_panel.structural_search_requested.connect(self._structural_search)
        self.pathway_panel.upstream_requested.connect(self._explore_upstream)
        # Thermodynamics is opt-in: hidden unless enabled in Settings ▸ Preferences.
        self.pathway_panel.set_mdf_enabled(self._mdf_enabled())
        self.pathway_panel.row_context_requested.connect(self._on_pathway_result_context)
        self.pathway_panel.manage_databases_requested.connect(self._manage_databases)
        # Enzyme (EC -> UniProt) lookup from the inspector
        self.info.enzymes_requested.connect(self._lookup_enzymes)
        # Categories
        self.categories_panel.new_requested.connect(self._new_category)
        self.categories_panel.delete_requested.connect(self._delete_category)
        self.categories_panel.add_selected_requested.connect(self._add_selected_to_category)
        self.categories_panel.remove_selected_requested.connect(self._remove_selected_from_category)
        self.categories_panel.isolate_requested.connect(self._isolate_category)
        self.categories_panel.analyze_requested.connect(self._analyze_category)

    def _new_menu(self, parent, title: str):
        """Add a menu (or submenu) that C++ owns — see ``_build_actions`` for why.

        Constructing the QMenu with its parent, rather than letting ``addMenu(title)``
        build one, is what moves ownership to C++. A Python-owned menu is destroyed as
        soon as the last wrapper referencing it is collected, and any later
        ``action.menu()`` makes a *fresh* Python-owned wrapper for the same object — so
        holding a reference is not enough on its own.
        """
        from PySide6.QtWidgets import QMenu

        menu = QMenu(title, parent)
        parent.addMenu(menu)
        self._menus.append(menu)
        return menu

    def _build_actions(self) -> None:
        # The menu bar and every menu are held on the window on purpose.
        #
        # ``QMainWindow.menuBar()`` creates the bar on first call and hands it to Python,
        # not to C++. Whenever the last Python reference goes away the underlying C++
        # QMenuBar is destroyed and takes every menu with it — so a throwaway
        # ``self.menuBar()`` anywhere in the codebase silently deletes the whole menu
        # bar, and the next call quietly builds an empty replacement. Until this was
        # pinned the menus survived only because nothing had triggered the collection
        # yet. The same applies to each QMenu returned by ``addMenu``.
        mb = self._menu_bar = self.menuBar()
        self._menus = []

        # File
        file_menu = self._new_menu(mb, "&File")
        self.act_open_model = QAction("Open Model…", self)
        self.act_open_model.setShortcut(QKeySequence.Open)
        self.act_open_model.triggered.connect(self.open_model)
        self.act_open_example = QAction("Open Example Model (E. coli core)", self)
        self.act_open_example.triggered.connect(self.open_example_model)
        self.act_open_project = QAction("Open Project…", self)
        self.act_open_project.triggered.connect(self.open_project)
        self.act_save_project = QAction("Save Project", self)
        self.act_save_project.setShortcut(QKeySequence.Save)
        self.act_save_project.triggered.connect(self.save_project)
        self.act_save_project_as = QAction("Save Project As…", self)
        self.act_save_project_as.triggered.connect(self.save_project_as)
        self.act_export_model = QAction("Export Model…", self)
        self.act_export_model.triggered.connect(self.export_model)
        self.act_quit = QAction("Quit", self)
        self.act_quit.triggered.connect(self.close)
        for a in (self.act_open_model, self.act_open_example, self.act_open_project,
                  self.act_save_project, self.act_save_project_as, self.act_export_model):
            file_menu.addAction(a)
        file_menu.addSeparator()
        file_menu.addAction(self.act_quit)

        # Edit
        edit_menu = self._new_menu(mb, "&Edit")
        self.act_undo = QAction("Undo", self)
        self.act_undo.setShortcut(QKeySequence.Undo)
        self.act_undo.triggered.connect(self.undo)
        self.act_redo = QAction("Redo", self)
        self.act_redo.setShortcut(QKeySequence.Redo)
        self.act_redo.triggered.connect(self.redo)
        self.act_add_rxn = QAction("Add Reaction…", self)
        self.act_add_rxn.triggered.connect(self.add_reaction)
        self.act_edit_rxn = QAction("Edit Selected Reaction…", self)
        self.act_edit_rxn.triggered.connect(self._edit_selected_reaction)
        self.act_remove_rxn = QAction("Remove Selected Reaction", self)
        self.act_remove_rxn.triggered.connect(self.remove_reaction)
        self.act_set_obj = QAction("Set Selected as Objective", self)
        self.act_set_obj.triggered.connect(self.set_objective)
        # Edit is Undo/Redo only. The reaction-editing actions still exist and are still
        # reachable — they live under Tools ▸ Edit Model, which is where a user looks for
        # "change the model" rather than in a menu shared with undo history.
        for a in (self.act_undo, self.act_redo):
            edit_menu.addAction(a)

        # Actions reused inside the Tools menu (built after View/Settings).
        self.act_media = QAction("Edit Growth Medium…", self)
        self.act_media.triggered.connect(self.edit_medium)
        self.act_build_community = QAction("Build Community Model…", self)
        self.act_build_community.triggered.connect(self.build_community)

        # View — which panels are on screen, and nothing else. Navigation ("Go to
        # Analysis") belongs to the tab strip, and panel *placement* is a setting, so
        # both have moved out.
        view_menu = self._new_menu(mb, "&View")
        view_menu.addAction(self.explorer_dock.toggleViewAction())
        view_menu.addAction(self.categories_dock.toggleViewAction())
        view_menu.addAction(self.info_dock.toggleViewAction())
        view_menu.addSeparator()
        self.act_model_info = QAction("Model info…", self)
        self.act_model_info.triggered.connect(self._show_model_info)
        view_menu.addAction(self.act_model_info)

        # Information-panel placement: a preference, so it lives under Settings.
        self.act_info_right = QAction("Right side", self)
        self.act_info_right.setCheckable(True)
        self.act_info_right.triggered.connect(lambda: self._set_info_position("right"))
        self.act_info_bottom = QAction("Bottom (under centre)", self)
        self.act_info_bottom.setCheckable(True)
        self.act_info_bottom.triggered.connect(lambda: self._set_info_position("bottom"))

        # Settings
        settings_menu = self._new_menu(mb, "&Settings")
        self.act_prefs = QAction("Preferences…", self)
        self.act_prefs.setToolTip("Every configurable option, grouped by category.")
        self.act_prefs.triggered.connect(self._open_preferences)
        settings_menu.addAction(self.act_prefs)
        settings_menu.addSeparator()
        info_pos = self._new_menu(settings_menu, "Information panel position")
        info_pos.addAction(self.act_info_right)
        info_pos.addAction(self.act_info_bottom)
        settings_menu.addSeparator()
        self.act_manage_data = QAction("Manage Data…", self)
        self.act_manage_data.setToolTip("View and clean the cached databases and molecule "
                                        "structure images stored on disk.")
        self.act_manage_data.triggered.connect(self._open_manage_data)
        settings_menu.addAction(self.act_manage_data)
        # The results folder is now a field in Preferences ▸ General, where it sits with
        # the other output options rather than as a lone menu entry.
        self.act_results_dir = QAction("Set results output folder…", self)
        self.act_results_dir.triggered.connect(self._set_results_dir)

        # Tools — a launcher for every utility (#T1)
        self._build_tools_menu(mb)

        # Help
        help_menu = self._new_menu(mb, "&Help")
        about = QAction("About", self)
        about.triggered.connect(self._about)
        help_menu.addAction(about)
        check_updates = QAction("Check for updates…", self)
        check_updates.triggered.connect(self._check_for_updates)
        help_menu.addAction(check_updates)

        self._build_toolbar()

    def _dynamic_environment(self) -> dict:
        """Baseline environment the regulatory sensors read during a dynamic run."""
        return {"light_uE": 300.0, "inorganic_c_mM": 50.0,
                "nitrogen_mM": 17.6, "iron_uM": 5.0}

    #: Environment key → (exchange reaction, mmol gDW⁻¹ h⁻¹ of uptake capacity per unit).
    #: The conversions are coarse: they translate a medium concentration into an uptake
    #: bound, which is a modelling choice rather than a measurement, and the panel says so.
    _DYNAMIC_MAPPING = {
        "light_uE": ("EX_photon_e", 0.18),      # a* = 0.05 m² gDW⁻¹ (see the model study)
        "inorganic_c_mM": ("EX_hco3_e", 0.5),
        "nitrogen_mM": ("EX_no3_e", 0.5),
        "iron_uM": ("EX_fe3_e", 0.1),
    }

    def _run_condition_scan(self, spec: dict) -> None:
        """Sweep one environmental variable and report how the network responds."""
        if self.project is None:
            QMessageBox.information(self, "No model", "Load a model first.")
            return
        import numpy as np
        from ..core import screening as scr

        from ..core import regulation as reg

        ruleset = self._active_ruleset(spec)

        low, high = float(spec["from"]), float(spec["to"])
        n = int(spec["points"])
        if spec.get("logarithmic"):
            low = max(low, 1e-4)
            grid = list(np.logspace(np.log10(low), np.log10(max(high, low * 10)), n))
        else:
            grid = list(np.linspace(low, high, n))

        model = self.project.model.copy()
        variable = spec["variable"]
        targets = spec.get("targets") or []

        # A scan over an exchange sweeps that exchange's uptake bound directly, and is
        # reported under the sensor it feeds so a rule set still reads it correctly.
        mapping = dict(self._DYNAMIC_MAPPING)
        if spec.get("kind") == "exchange":
            if not model.reactions.has_id(variable):
                QMessageBox.warning(self, "Not in this model",
                                    f"This model has no {variable}.")
                return
            sensor = reg.sensor_for_exchange(variable)
            mapping[sensor] = (variable, 1.0)
            variable = sensor

        def work():
            return scr.scan(model, ruleset,
                            variables={variable: grid},
                            base_environment=self._dynamic_environment(),
                            mapping=mapping,
                            targets=targets)

        self.dynamic_panel.set_busy(True)
        try:
            ok, result = run_busy(self, f"Scanning {spec.get('label') or variable}…",
                                  work, title="Condition scan", cancelable=True)
        finally:
            self.dynamic_panel.set_busy(False)
        if not ok:
            if not was_cancelled(result):
                QMessageBox.warning(self, "Scan failed", str(result))
            return

        import pandas as pd
        frame = pd.DataFrame(result.table())
        breaks = result.transitions()
        commentary = [f"{len(grid)} points scanned across {variable}."]
        if breaks:
            commentary.append("Breakpoints — where the response changes sharply:")
            for b in breaks:
                a = list(b["from"].values())[0]
                z = list(b["to"].values())[0]
                commentary.append(
                    f"  • between {a:.4g} and {z:.4g}: growth "
                    f"{b['growth'][0]:.5f} → {b['growth'][1]:.5f} "
                    f"({b['relative_change'] * 100:.0f}% change)")
        else:
            commentary.append("No sharp breakpoint — the response is smooth across this "
                              "range.")
        warning = self._regulation_warning(result.points)
        self._dynamic_frame = frame
        self._dynamic_title = f"Condition scan — {spec.get('label') or variable}"
        self.dynamic_panel.show_table(frame, self._dynamic_title,
                                      "\n".join(commentary), warning, kind="scan")
        self.tabs.setCurrentWidget(self.dynamic_panel)

    def _active_ruleset(self, spec: dict):
        """The rule set a dynamic run should use, honouring the panel's checkbox.

        Regulation stays opt-in: with it off, an empty rule set makes the run behave
        exactly as the unregulated model, which is what makes the feature safe to ship.
        """
        from ..core import preferences as prefs
        from ..core import regulation as reg

        if not spec.get("regulation", prefs.get(prefs.REGULATION_ENABLED)):
            return reg.RuleSet()
        ruleset, _path = self._regulation_ruleset()
        return ruleset

    def _run_timecourse(self, spec: dict) -> None:
        """Follow a batch culture as its chosen medium components deplete."""
        if self.project is None:
            QMessageBox.information(self, "No model", "Load a model first.")
            return
        from ..core import screening as scr

        ruleset = self._active_ruleset(spec)
        model = self.project.model.copy()

        chosen = spec.get("substrates") or []
        if not chosen:
            QMessageBox.information(
                self, "Nothing to follow",
                "Add at least one medium component to the table. A time course is the "
                "story of something running out — without a component to track there is "
                "nothing to integrate.")
            return

        substrates, missing = [], []
        for item in chosen:
            exchange = item["exchange"]
            if not model.reactions.has_id(exchange):
                missing.append(exchange)
                continue
            substrates.append(scr.substrate_from_exchange(
                model, exchange, item["initial_mM"],
                sensor=item.get("sensor", ""),
                max_uptake=item.get("max_uptake", 0.0),
                buffered=bool(item.get("buffered"))))
        if missing:
            QMessageBox.warning(
                self, "Components not in model",
                "This model has no " + ", ".join(missing) +
                ".\n\nThey were skipped; the run continues with the rest.")
        if not substrates:
            return

        environment = self._dynamic_environment()
        environment["light_uE"] = float(spec["light_uE"])
        targets = spec.get("targets") or []

        def work():
            return scr.timecourse(
                model, ruleset, substrates=substrates,
                base_environment=environment,
                initial_biomass=float(spec["biomass"]),
                duration_h=float(spec["duration_h"]),
                step_h=float(spec["step_h"]),
                targets=targets)

        label = ", ".join(s.column() for s in substrates[:3])
        self.dynamic_panel.set_busy(True)
        try:
            ok, result = run_busy(self, f"Simulating {label}…", work,
                                  title="Time course", cancelable=True)
        finally:
            self.dynamic_panel.set_busy(False)
        if not ok:
            if not was_cancelled(result):
                QMessageBox.warning(self, "Time course failed", str(result))
            return

        import pandas as pd
        frame = pd.DataFrame(result.table())
        commentary = list(result.notes)
        depleted = [s.column() for s in substrates
                    if not s.buffered
                    and result.points and
                    result.points[-1].concentrations.get(s.column(), 1.0) <= 1e-6]
        if depleted:
            commentary.append("Exhausted by the end of the run: " + ", ".join(depleted)
                              + ".")
        changes = result.phase_changes()
        if changes:
            commentary.append("Regulatory phase changes:")
            for c in changes:
                gained = ", ".join(c["gained"]) or "—"
                commentary.append(
                    f"  • t = {c['time_h']:.0f} h: {gained} "
                    f"(growth {c['growth_before']:.4f} → {c['growth_after']:.4f} h⁻¹)")
        elif len(ruleset):
            commentary.append("No regulatory transition occurred during this run.")
        else:
            commentary.append("Run without regulation — no rules were applied.")

        self._dynamic_frame = frame
        self._dynamic_title = f"Time course — {label}"
        self.dynamic_panel.show_table(frame, self._dynamic_title,
                                      "\n".join(commentary), kind="timecourse")
        self.tabs.setCurrentWidget(self.dynamic_panel)

    def _plot_dynamic_results(self, names=None) -> None:
        """Open the axis chooser over the stored runs, ready to overlay them."""
        runs = self.dynamic_panel.runs_for_plot()
        if names:
            wanted = set(names)
            runs = [(n, f) for n, f in runs if n in wanted] or runs
        runs = [(n, f) for n, f in runs if f is not None and not f.empty]
        if not runs:
            QMessageBox.information(self, "Nothing to plot",
                                    "Run a scan or a time course first.")
            return
        from .dialogs.timecourse_plot_dialog import TimecoursePlotDialog

        current = self.dynamic_panel.current_run()
        title = current or getattr(self, "_dynamic_title", "Dynamic run")
        # Start on the selected run; the dialog offers the rest for overlay.
        ordered = sorted(runs, key=lambda item: item[0] != current)
        TimecoursePlotDialog(self, ordered[0][1], title, runs=ordered).exec()

    def _show_dynamic_run(self, name: str) -> None:
        """Open one stored run's table in its own window."""
        record = self.dynamic_panel.run_record(name)
        if record is None:
            return
        from .dialogs.run_results_dialog import RunResultsDialog

        dialog = RunResultsDialog(self, name, record)
        dialog.plot_requested.connect(lambda n: self._plot_dynamic_results([n]))
        dialog.show()

    @staticmethod
    def _regulation_warning(points) -> str:
        """Flag when a scan's conclusion rests on rules with assumed thresholds."""
        from ..core import regulation as reg

        worst = reg.MEASURED
        order = {reg.MEASURED: 0, reg.INFERRED: 1, reg.ASSUMED: 2}
        for p in points:
            if order.get(p.weakest_confidence, 0) > order.get(worst, 0):
                worst = p.weakest_confidence
        if worst == reg.ASSUMED:
            return ("<span style='color:#c5221f'>⚠ Some active regulatory rules use "
                    "<b>assumed</b> thresholds — read these numbers as indicative.</span>")
        if worst == reg.INFERRED:
            return ("<span style='color:#b06000'>Some active rules use <b>inferred</b> "
                    "thresholds.</span>")
        return ""

    def _focus_dynamic(self, index: int) -> None:
        self.tabs.setCurrentWidget(self.dynamic_panel)
        for child in self.dynamic_panel.findChildren(QTabWidget):
            child.setCurrentIndex(index)
            break

    def _regulation_ruleset(self):
        """The rule set the user made active, or an empty one.

        Empty is the honest default. Nothing is bundled as a fallback: a rule set encodes
        one organism's regulation, so quietly supplying somebody else's would produce
        confident answers about the wrong biology.
        """
        from ..core import rule_library as lib

        return lib.active()

    def _open_regulation(self) -> None:
        """Tools ▸ Regulation — the rule-set library."""
        from .dialogs.regulation_dialog import RegulationDialog

        model = self.project.model if self.project is not None else None
        RegulationDialog(self, model).exec()
        self._sync_regulation_status()

    def _reload_ruleset(self) -> None:
        ruleset, path = self._regulation_ruleset()
        if not path:
            self.status_label.setText(
                "No regulatory rule set is active — simulations run unregulated.")
        else:
            self.status_label.setText(
                f"Reloaded {len(ruleset)} regulatory rule(s) from "
                f"{os.path.basename(path)}.")
        self._sync_regulation_status()

    def _sync_regulation_status(self) -> None:
        """Tell the Dynamic Analysis tab which rule set is in force, if any."""
        from ..core import preferences as prefs

        panel = getattr(self, "dynamic_panel", None)
        if panel is None:
            return
        ruleset, path = self._regulation_ruleset()
        if not path:
            panel.set_regulation_status(
                "No rule set active — runs are unregulated. "
                "Use <i>Regulatory rules…</i> to load or write one.", enabled=False)
            panel.use_regulation.setEnabled(False)
            return
        panel.use_regulation.setEnabled(True)
        organism = getattr(ruleset, "organism", "") or "organism not recorded"
        panel.set_regulation_status(
            f"Active: <b>{ruleset.name or os.path.basename(path)}</b> — "
            f"{len(ruleset.enabled())} rule(s), {organism}.",
            enabled=bool(prefs.get(prefs.REGULATION_ENABLED)))

    def _set_info_position(self, where: str) -> None:
        """Place the Information dock, and remember the choice."""
        from ..core import preferences as prefs

        area = Qt.RightDockWidgetArea if where == "right" else Qt.BottomDockWidgetArea
        with self.hold_window_state():
            self.addDockWidget(area, self.info_dock)
            self.info_dock.show()
        prefs.set(prefs.INFO_PANEL_POSITION, where)
        self.act_info_right.setChecked(where == "right")
        self.act_info_bottom.setChecked(where == "bottom")

    def _toolbar_action(self, action_id: str):
        """Map a stored shortcut id onto a QAction, or None if it is not available.

        Ids are stored rather than actions so a saved bar survives the actions being
        rebuilt, and an id that no longer exists is skipped instead of crashing.
        """
        simple = {
            "open_model": self.act_open_model,
            "save_project": self.act_save_project,
            "add_reaction": self.act_add_rxn,
            "undo": self.act_undo,
            "redo": self.act_redo,
            "edit_medium": self.act_media,
            "preferences": self.act_prefs,
        }
        if action_id in simple:
            return simple[action_id]

        if action_id == "growth_settings":
            act = QAction("Growth Settings", self)
            act.setToolTip("Light, medium and uptake settings for the loaded model.")
            act.triggered.connect(self._open_growth_settings)
            return act
        if action_id == "network_visualization":
            return None            # handled separately: it is a button with a menu
        if action_id == "pathway_design":
            act = QAction("Pathway Design", self)
            act.triggered.connect(lambda: self.tabs.setCurrentWidget(self.pathway_panel))
            return act
        if action_id == "dynamic_analysis":
            act = QAction("Dynamic Analysis", self)
            act.triggered.connect(lambda: self.tabs.setCurrentWidget(self.dynamic_panel))
            return act
        if action_id == "manage_databases":
            act = QAction("Manage Databases", self)
            act.triggered.connect(self._manage_databases)
            return act
        if action_id in ("run_fba", "run_pfba"):
            label = "Run FBA" if action_id == "run_fba" else "Run pFBA"
            act = QAction(label, self)
            aid = "fba" if action_id == "run_fba" else "pfba"
            act.triggered.connect(lambda _=False, a=aid: self._run_analysis(a))
            return act
        return None

    def _build_toolbar(self) -> None:
        """The quick-access bar, built from the user's chosen shortcut ids."""
        from ..core import preferences as prefs

        if getattr(self, "_toolbar", None) is None:
            self._toolbar = self.addToolBar("Quick access")
            self._toolbar.setMovable(False)
        self._toolbar.clear()
        self._toolbar_extras = []          # keep dynamically created actions alive

        chosen = prefs.get(prefs.TOOLBAR_ACTIONS) or list(prefs.TOOLBAR_DEFAULT)
        for action_id in chosen:
            if action_id == "network_visualization":
                self._toolbar.addWidget(self._network_visualization_button())
                continue
            action = self._toolbar_action(action_id)
            if action is not None:
                self._toolbar.addAction(action)
                self._toolbar_extras.append(action)

    def _network_visualization_button(self):
        """One button for the three network views, instead of three main tabs.

        Escher and Strategy were full tabs competing for the tab strip while being used
        occasionally. Collapsing them into a menu button frees the strip for Dynamic
        Analysis without hiding anything.
        """
        from PySide6.QtWidgets import QMenu, QToolButton

        button = QToolButton()
        button.setText("Network Visualization")
        button.setPopupMode(QToolButton.InstantPopup)
        button.setToolTip("Open a network view: interactive Escher map, strategy "
                          "comparison, or the static network map.")
        menu = QMenu(button)
        menu.addAction("Escher Visualizer", self._open_escher)
        menu.addAction("Strategy Visualizer",
                       lambda: self._open_floating(self.strategy_explorer,
                                                   "Strategy Visualizer"))
        menu.addSeparator()
        menu.addAction("Network Map",
                       lambda: self.tabs.setCurrentWidget(self.network_view))
        button.setMenu(menu)
        button.setEnabled(self.project is not None)
        self._network_menu = menu           # keep both alive (see _build_actions)
        self._network_button = button
        return button

    def _open_floating(self, widget, title: str) -> None:
        """Show a panel that no longer has a tab, in its own window."""
        from PySide6.QtWidgets import QDialog, QVBoxLayout

        holder = getattr(self, "_floating_windows", None)
        if holder is None:
            holder = self._floating_windows = {}
        existing = holder.get(title)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return
        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog.resize(1000, 700)
        layout = QVBoxLayout(dialog)
        layout.addWidget(widget)
        from .widgets.dialog_util import clamp_to_screen
        clamp_to_screen(dialog)
        holder[title] = dialog
        dialog.show()

    def _open_growth_settings(self) -> None:
        """Growth settings, now a toolbar button rather than a permanent tab."""
        self._open_floating(self.growth_panel, "Growth Settings")

    def _build_tools_menu(self, mb) -> None:
        """Build the top-level Tools menu: every utility as menus/submenus (#T1)."""
        from .panels.analysis_panel import (SIMULATION, STRAIN_DESIGN, PHENOTYPE,
                                            ESSENTIALITY, OMICS, PATHWAYS, CURATION)
        tools = self._new_menu(mb, "&Tools")

        def _go(widget):
            return lambda: self.tabs.setCurrentWidget(widget)

        def _analysis(aid):
            def _do():
                self.tabs.setCurrentWidget(self.analysis_panel)
                self.analysis_panel.run_requested.emit(aid)
            return _do

        # Model-dependent submenus, disabled together when nothing is loaded.
        self._tools_model_menus = []

        edit = self._new_menu(tools, "Edit Model")
        edit.addAction("Add reaction…", self.add_reaction)
        edit.addAction(self.act_media)                       # Edit growth medium…
        edit.addAction("Growth settings…", self._open_growth_settings)
        edit.addAction("Create category…", self._new_category)
        edit.addSeparator()
        edit.addAction("Detect subsystems…", self._detect_subsystems)
        self._tools_model_menus.append(edit)

        analysis = self._new_menu(tools, "Analysis")
        for label, group in (("Simulation", SIMULATION), ("Strain design", STRAIN_DESIGN),
                             ("Phenotype & essentiality", PHENOTYPE + ESSENTIALITY),
                             ("Omics & energy", OMICS), ("Pathways & community", PATHWAYS),
                             ("Model quality & curation", CURATION)):
            sub = self._new_menu(analysis, label)
            for aid, lbl, tip in group:
                act = sub.addAction(lbl.replace("\n", " "), _analysis(aid))
                act.setToolTip(tip)
        self._tools_model_menus.append(analysis)

        # Pathway Design stays enabled without a model: reaction databases can be
        # downloaded and merged before there is anything to apply them to.
        pw = self._new_menu(tools, "Pathway Design")
        pw.addAction("Open Pathway Design", _go(self.pathway_panel))
        pw.addAction("Manage reaction databases…", self._manage_databases)
        pw.addAction("Reconcile metabolite identifiers…", self._reconcile_identifiers)

        self.act_community = tools.addAction(
            "Consortia modelling — Build Community Model…", self.build_community)

        viz = self._new_menu(tools, "Visualization")
        viz.addAction("Network Map", _go(self.network_view))
        viz.addAction("Strategy Visualizer",
                      lambda: self._open_floating(self.strategy_explorer,
                                                  "Strategy Visualizer"))
        viz.addAction("Escher Visualizer", self._open_escher)   # interactive Escher maps (#T6)
        self._tools_model_menus.append(viz)

        dyn = self._new_menu(tools, "Dynamic Analysis")
        dyn.addAction("Open Dynamic Analysis", _go(self.dynamic_panel))
        dyn.addAction("Condition scan", lambda: self._focus_dynamic(0))
        dyn.addAction("Time course (batch)", lambda: self._focus_dynamic(1))
        self._tools_model_menus.append(dyn)

        # Regulation gets its own submenu rather than being buried in Edit Model: it
        # changes how every simulation behaves, so it deserves to be found.
        regulation = self._new_menu(tools, "Regulation")
        regulation.addAction("Regulatory model…", self._open_regulation)
        regulation.addAction("Reload rule set", self._reload_ruleset)

    def _build_statusbar(self) -> None:
        self.objective_label = QLabel("")
        # Job-aware status: specific description + estimated-progress bar + expand
        # arrow that opens a per-process dialog (handles several concurrent jobs).
        self.status_label = JobStatusBar(self.jobs)
        self.statusBar().addWidget(self.status_label, 1)
        self.statusBar().addPermanentWidget(self.objective_label)

    # ----- model / project lifecycle -----------------------------------
    def open_model(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open metabolic model", "", _MODEL_FILTER)
        if not path:
            return
        self._load_project(lambda: Project.from_model_file(path), os.path.basename(path))

    def open_example_model(self) -> None:
        from .. import resources

        path = resources.example_model_path()
        self._load_project(lambda: Project.from_model_file(path), "E. coli core")

    def open_project(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open project", "", _PROJECT_FILTER)
        if not path:
            return
        self._load_project(lambda: Project.load(path), os.path.basename(path))

    def _load_project(self, factory: Callable[[], Project], label: str) -> None:
        # Parse the model on a worker thread, then populate all the panels
        # (explorer tables, targets, etc.) *while the loading dialog is still
        # shown* — so it stays up until the app is genuinely ready to use rather
        # than vanishing and leaving the window frozen mid-refresh.
        def _finish(project):
            self.project = project
            self._last_result = None
            self.pathway_panel.clear_results()   # a new model invalidates old pathways
            self._refresh_all()

        ok, result = run_busy(self, f"Loading {label}…", factory, title="Loading model",
                              after=_finish, after_message=f"Preparing {label}…",
                              cancelable=True)
        if not ok:
            if not was_cancelled(result):
                QMessageBox.critical(self, "Could not open", str(result))
            return
        self.status_label.setText(f"Loaded {label}.")
        self._show_model_summary_on_load()

    def save_project(self) -> None:
        if self.project is None:
            return
        if not self.project.project_path:
            self.save_project_as()
            return
        try:
            self.project.save()
        except ProjectError as exc:
            QMessageBox.critical(self, "Could not save project", str(exc))
            return
        self._update_title()
        self.status_label.setText("Project saved.")

    def save_project_as(self) -> None:
        if self.project is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save project", "", _PROJECT_FILTER)
        if not path:
            return
        try:
            self.project.save(path)
        except ProjectError as exc:
            QMessageBox.critical(self, "Could not save project", str(exc))
            return
        self._update_title()
        self.status_label.setText("Project saved.")

    def export_model(self) -> None:
        if self.project is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export model", "", _MODEL_FILTER)
        if not path:
            return
        try:
            io_models.save_model(self.project.model, path)
        except io_models.ModelSaveError as exc:
            QMessageBox.critical(self, "Could not export model", str(exc))
            return
        self.status_label.setText(f"Exported model to {os.path.basename(path)}.")

    # ----- editing ------------------------------------------------------
    def add_reaction(self) -> None:
        if self.project is None:
            return
        dlg = ReactionBuilderDialog(self.project.model, self,
                                    database=getattr(self, "_pathway_db", None))
        if dlg.exec() != ReactionBuilderDialog.Accepted:
            return
        v = dlg.values()
        ec = v.get("ec_number", "")
        met_info = v.get("metabolite_info", {})

        def _apply(m):
            rxn = editing.add_reaction(
                m, v["reaction_id"], v["reaction_string"], name=v["name"],
                lower_bound=v["lower_bound"], upper_bound=v["upper_bound"],
                gene_reaction_rule=v["gene_reaction_rule"], subsystem=v["subsystem"])
            if ec:
                rxn.annotation["ec-code"] = [e.strip() for e in ec.split(";") if e.strip()]
            # Give newly-created metabolites proper names + cross-references so
            # they display correctly (and match across databases later).
            for met in rxn.metabolites:
                info = met_info.get(met.id)
                if info:
                    if info.get("name") and not met.name:
                        met.name = info["name"]
                    for key, val in (info.get("annotation") or {}).items():
                        met.annotation.setdefault(key, val)
                    if not met.compartment and met.id.rsplit("_", 1)[-1]:
                        met.compartment = met.id.rsplit("_", 1)[-1]
            return rxn

        ok, result = run_busy(self, f"Adding reaction {v['reaction_id']}…",
                              lambda: self.project.apply_edit(_apply), title="Adding reaction")
        if not ok:
            QMessageBox.warning(self, "Could not add reaction", str(result))
            return
        self._refresh_after_edit(f"Added reaction {v['reaction_id']}.")

    def _edit_selected_reaction(self) -> None:
        rxn = self.explorer.selected_reaction()
        if rxn is None:
            QMessageBox.information(self, "No selection", "Select a reaction first.")
            return
        self._edit_reaction(rxn)

    def _edit_reaction(self, rxn) -> None:
        if self.project is None:
            return
        dlg = ReactionDialog(self, reaction=rxn)
        if dlg.exec() != ReactionDialog.Accepted:
            return
        v = dlg.values()
        rid = rxn.id

        def _apply(m):
            editing.set_bounds(m, rid, v["lower_bound"], v["upper_bound"])
            editing.set_gene_reaction_rule(m, rid, v["gene_reaction_rule"])
            r = m.reactions.get_by_id(rid)
            r.name = v["name"]
            r.subsystem = v["subsystem"]
            if v.get("balance_changes"):     # staged H+/H2O balancing (fix 12)
                from ..core import balancing
                balancing.apply_changes(r, v["balance_changes"])

        try:
            self.project.apply_edit(_apply)
        except editing.EditError as exc:
            QMessageBox.warning(self, "Could not edit reaction", str(exc))
            return
        self._refresh_after_edit(f"Updated reaction {rid}.")

    def remove_reaction(self) -> None:
        if self.project is None:
            return
        rxn = self.explorer.selected_reaction()
        if rxn is None:
            QMessageBox.information(self, "No selection", "Select a reaction to remove.")
            return
        rid = rxn.id
        if QMessageBox.question(self, "Remove reaction", f"Remove reaction '{rid}'?") != QMessageBox.Yes:
            return
        self.project.apply_edit(lambda m: editing.remove_reaction(m, rid))
        self._refresh_after_edit(f"Removed reaction {rid}.")

    def set_objective(self) -> None:
        if self.project is None:
            return
        rxn = self.explorer.selected_reaction()
        if rxn is None:
            QMessageBox.information(self, "No selection", "Select a reaction to set as objective.")
            return
        rid = rxn.id
        self.project.apply_edit(lambda m: editing.set_objective(m, rid, "max"))
        self._refresh_after_edit(f"Objective set to maximize {rid}.")

    def edit_medium(self) -> None:
        if self.project is None:
            return
        dlg = MediaDialog(self.project.model, self)
        if dlg.exec() != MediaDialog.Accepted:
            return
        new_medium = dlg.medium()
        aerobic = dlg.is_aerobic()

        def _apply(m):
            media.set_medium(m, new_medium)
            media.set_aerobic(m, aerobic)

        self.project.apply_edit(_apply)
        self._refresh_after_edit("Updated growth medium.")

    def undo(self) -> None:
        if self.project and self.project.can_undo:
            self.project.undo()
            self._last_result = None
            self._refresh_all()

    def redo(self) -> None:
        if self.project and self.project.can_redo:
            self.project.redo()
            self._last_result = None
            self._refresh_all()

    # ----- categories ---------------------------------------------------
    def _detect_subsystems(self) -> None:
        """Tools ▸ Edit Model ▸ Detect subsystems.

        Three steps on purpose: ask what to detect, detect, then review before anything
        is written. Detection is inference — writing it straight onto the model would
        replace an honest blank with a confident label nobody checked, and every feature
        that groups by subsystem would inherit the mistake without a trace.
        """
        if self.project is None:
            QMessageBox.information(self, "No model", "Load a model first.")
            return
        from ..core import subsystems as subs
        from .dialogs.subsystem_dialog import SubsystemDialog, SubsystemOptionsDialog

        model = self.project.model
        options = SubsystemOptionsDialog(self, model)
        if not options.exec():
            return
        settings = options.values()

        def work():
            return subs.detect(model, overwrite=settings["overwrite"],
                               minimum_evidence=settings["minimum_evidence"])

        ok, report = run_busy(self, "Matching reactions to central pathways…", work,
                              title="Detect subsystems")
        if not ok:
            if not was_cancelled(report):
                QMessageBox.warning(self, "Detection failed", str(report))
            return
        if not report.assignments:
            QMessageBox.information(
                self, "Nothing detected",
                "No reaction matched a central pathway at this evidence level.\n\n"
                "Try again accepting weaker evidence, or check that the model carries "
                "recognisable reaction ids or EC annotations.")
            return

        review = SubsystemDialog(self, model, report)
        if not review.exec():
            return
        assignments = review.assignments()

        def _apply(m):
            return subs.apply(m, assignments)

        changed = self.project.apply_edit(_apply)
        self._refresh_all()
        self.status_label.setText(
            f"Subsystems written to {changed} reaction(s) — undo is available.")
        QMessageBox.information(
            self, "Subsystems applied",
            f"{changed} reaction(s) annotated across "
            f"{len(set(assignments.values()))} subsystem(s).\n\n"
            "The Explorer, Escher focus and pathway grouping now use them. Use Edit ▸ "
            "Undo to revert.")

    def _new_category(self) -> None:
        if self.project is None:
            return
        name, ok = QInputDialog.getText(self, "New category", "Category name:")
        if not ok or not name.strip():
            return
        try:
            self.project.categories.create(name.strip())
        except ValueError as exc:
            QMessageBox.warning(self, "Category", str(exc))
            return
        self._refresh_categories()

    def _delete_category(self, name: str) -> None:
        if self.project is None:
            return
        self.project.categories.delete(name)
        self._refresh_categories()

    def _add_selected_to_category(self, name: str) -> None:
        ids = self.explorer.selected_reaction_ids()
        if not ids:
            QMessageBox.information(self, "No selection",
                                    "Select one or more reactions in the Explorer first.")
            return
        self._add_ids_to_category(name, ids)

    def _add_ids_to_category(self, name: str, ids) -> None:
        if self.project is None or not self.project.categories.has(name):
            return
        self.project.categories.get(name).add(ids)
        self.project._dirty = True
        self._refresh_categories()
        self.status_label.setText(f"Added {len(ids)} reaction(s) to '{name}'.")

    def _add_reactions_to_category(self, ids: list) -> None:
        """From the explorer context menu: choose/create a category, then add `ids`."""
        if self.project is None or not ids:
            return
        names = self.project.categories.names()
        new_label = "New category…"
        choice = new_label
        if names:
            item, ok = QInputDialog.getItem(
                self, "Add to category", "Choose a category:", names + [new_label], 0, False)
            if not ok:
                return
            choice = item
        if choice == new_label:
            name, ok = QInputDialog.getText(self, "New category", "Category name:")
            if not ok or not name.strip():
                return
            name = name.strip()
            if not self.project.categories.has(name):
                self.project.categories.create(name)
        else:
            name = choice
        self._add_ids_to_category(name, ids)

    def _remove_selected_from_category(self, name: str) -> None:
        ids = self.explorer.selected_reaction_ids()
        if self.project and ids and self.project.categories.has(name):
            self.project.categories.get(name).remove(ids)
            self.project._dirty = True
            self._refresh_categories()

    def _isolate_category(self, name: str) -> None:
        if self.project is None or not self.project.categories.has(name):
            return
        ids = self.project.categories.get(name).reaction_ids
        if not ids:
            QMessageBox.information(self, "Empty category", "This category has no reactions yet.")
            return
        self.tabs.setCurrentWidget(self.network_view)
        self.network_view.focus_category(name)

    def _analyze_category(self, name: str) -> None:
        idx = self.analysis_panel.scope_combo.findData(name)
        if idx >= 0:
            self.analysis_panel.scope_combo.setCurrentIndex(idx)
        self.tabs.setCurrentWidget(self.analysis_panel)

    def _category_colors(self) -> dict:
        colors = {}
        if self.project is not None:
            for cat in self.project.categories.all():
                for rid in cat.reaction_ids:
                    colors[rid] = cat.color
        return colors

    def _refresh_categories(self) -> None:
        if self.project is None:
            return
        self.categories_panel.refresh(self.project.categories)
        self.analysis_panel.set_categories(self.project.categories.names())
        self.network_view.set_categories(self.project.categories.all())
        self.escher_explorer.set_categories(
            {c.name: sorted(c.reaction_ids) for c in self.project.categories.all()})

    # ----- map navigation ----------------------------------------------
    def _show_in_map(self, reaction_id: str, steps: int) -> None:
        self.tabs.setCurrentWidget(self.network_view)
        self.network_view.focus_on(reaction_id, steps)

    def _on_node_context(self, label: str, kind: str, global_pos) -> None:
        """Right-click menu on a map node (reaction or metabolite)."""
        if self.project is None:
            return
        model = self.project.model
        menu = QMenu(self)
        show_menu = menu.addMenu("Show in map")
        for steps in (1, 2, 3):
            act = show_menu.addAction(f"{steps} step{'s' if steps > 1 else ''} around")
            act.triggered.connect(lambda _=False, s=steps, n=label: self.network_view.focus_on(n, s))

        if kind == "reaction" and model.reactions.has_id(label):
            rxn = model.reactions.get_by_id(label)
            menu.addAction("Edit reaction…", lambda: self._edit_reaction(rxn))
            menu.addAction("Set as objective", lambda: self._set_objective_id(label))
            menu.addSeparator()
            menu.addAction("Add to category…", lambda: self._add_reactions_to_category([label]))
        elif kind == "metabolite" and model.metabolites.has_id(label):
            met = model.metabolites.get_by_id(label)
            menu.addAction("Show details", lambda: self._show_metabolite_details(met))
            menu.addAction("Edit metabolite…", lambda: self._edit_metabolite(met))
            connected = [r.id for r in met.reactions]
            menu.addAction(f"Add {len(connected)} connected reaction(s) to category…",
                           lambda: self._add_reactions_to_category(connected))
        menu.exec(global_pos)

    def _set_objective_id(self, reaction_id: str) -> None:
        if self.project is None:
            return
        self.project.apply_edit(lambda m: editing.set_objective(m, reaction_id, "max"))
        self._refresh_after_edit(f"Objective set to maximize {reaction_id}.")

    def _apply_objective(self, terms: dict, direction: str) -> None:
        if self.project is None:
            return
        try:
            self.project.apply_edit(
                lambda m: editing.set_weighted_objective(m, terms, direction))
        except editing.EditError as exc:
            QMessageBox.warning(self, "Objective", str(exc))
            return
        pretty = " + ".join(f"{w:g}·{rid}" for rid, w in terms.items())
        self._refresh_after_edit(f"Objective set to {direction} {pretty}.")

    def _on_network_node(self, label: str) -> None:
        if self.project is None:
            return
        model = self.project.model
        if model.reactions.has_id(label):
            self.explorer.select_reaction_by_id(label)
            flux = None
            if self._last_result is not None and self._last_result.is_optimal:
                flux = self._last_result.fluxes.get(label)
            self.info.show_reaction(model.reactions.get_by_id(label), flux)
        elif model.metabolites.has_id(label):
            self.info.show_metabolite(model.metabolites.get_by_id(label))

    # ----- analyses (always on a thread-safe copy) ---------------------
    def _scoped_model(self):
        """Return a private model copy for analysis, honoring the category scope.

        Returns (model_copy, scope_label). For a category scope, builds an isolated
        sub-model with free boundary exchanges for cut metabolites.
        """
        cat_name = self.analysis_panel.current_category()
        if cat_name and self.project.categories.has(cat_name):
            ids = self.project.categories.get(cat_name).reaction_ids
            return cat_core.build_subset_model(self.project.model, ids), f"category '{cat_name}'"
        return self.project.model.copy(), "whole model"

    _ANALYSIS_TITLES = {
        "fva": "Flux Variability Analysis settings",
        "production_envelope": "Production envelope settings",
        "robustness": "Robustness analysis settings",
        "phase_plane": "Phenotypic phase plane settings",
        "quality_report": "Quality report settings",
        "mutant": "Mutant simulation settings",
        "knockout": "Knockout strain design settings",
        "overproduction": "Metabolite overproduction settings",
        "fseof": "FSEOF settings",
        "gimme": "GIMME settings",
        "atpm_sensitivity": "ATP maintenance sensitivity settings",
        "efm": "Elementary flux modes settings",
        "gapfill_growth": "Gap-fill settings",
        "gapfill_metabolite": "Gap-fill settings",
    }

    def _configure(self, analysis_id: str):
        """Open the per-analysis settings dialog; return values dict, or None if cancelled."""
        from .dialogs.analysis_config import (
            AnalysisConfigDialog,
            description_for,
            params_for,
        )

        sel = self.explorer.selected_reaction()
        context = {"selected_reaction": sel.id if sel is not None else None}
        expr = self.project.datasets.get("expression")
        if expr:
            import numpy as np
            context["expression_median"] = float(np.median(list(expr.values())))
        params = params_for(analysis_id, self.project.model, context)
        if not params:
            return {}
        dlg = AnalysisConfigDialog(
            self._ANALYSIS_TITLES.get(analysis_id, "Analysis settings"),
            params, self.project.model, self, description=description_for(analysis_id))
        if dlg.exec() != AnalysisConfigDialog.Accepted:
            return None
        return dlg.values()

    def _compute_exclusions(self, cfg: dict) -> None:
        """Record which transport/exchange reactions to hide from result tables."""
        from ..core.network_graph import reaction_type
        ids = set()
        if self.project is not None and (cfg.get("exclude_transport") or cfg.get("exclude_exchange")):
            for r in self.project.model.reactions:
                t = reaction_type(r)
                if (t == "transport" and cfg.get("exclude_transport")) or \
                   (t == "exchange" and cfg.get("exclude_exchange")):
                    ids.add(r.id)
        self._excluded_rxn_ids = ids

    def _run_analysis(self, analysis_id: str) -> None:
        if self.project is None:
            QMessageBox.information(self, "No model", "Open a model first.")
            return
        if analysis_id == "prepare_omics":
            self._prepare_omics_dataset()
            return
        if analysis_id == "load_expression":
            self._load_expression()
            return
        if analysis_id in ("eflux", "gimme") and not self.project.datasets.get("expression"):
            QMessageBox.information(
                self, "No expression data",
                "Load a gene-expression table first (Omics & energy ▸ Load expression data).")
            return

        cfg = self._configure(analysis_id)
        if cfg is None:  # user cancelled the settings dialog
            return
        self._compute_exclusions(cfg)

        if analysis_id == "efm":
            self._run_efm(cfg)
            return
        if analysis_id == "community_growth":
            self._run_community_growth(cfg)
            return
        try:
            model, scope = self._scoped_model()
        except ValueError as exc:
            QMessageBox.warning(self, "Scope error", str(exc))
            return

        handlers = {
            "fba": self._run_fba_like,
            "pfba": self._run_fba_like,
            "shadow_prices": self._run_fba_like,
            "fva": self._run_fva,
            "mutant": self._run_mutant,
            "single_reaction_deletion": self._run_single_reaction_deletion,
            "single_gene_deletion": self._run_single_gene_deletion,
            "production_envelope": self._run_production_envelope,
            "robustness": self._run_robustness,
            "phase_plane": self._run_phase_plane,
            "quality_report": self._run_quality_report,
            "blocked_reactions": self._run_blocked_reactions,
            "gapfill_growth": self._run_gapfill_growth,
            "gapfill_metabolite": self._run_gapfill_metabolite,
            "knockout": self._run_knockout_design,
            "overproduction": self._run_overproduction,
            "fseof": self._run_fseof,
            "eflux": self._run_eflux,
            "gimme": self._run_gimme,
            "atpm_sensitivity": self._run_atpm,
            "flux_sampling": self._run_flux_sampling,
        }
        handler = handlers.get(analysis_id)
        if handler is not None:
            handler(analysis_id, model, scope, cfg)

    def _run_fba_like(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        method = {"fba": "FBA", "pfba": "pFBA", "shadow_prices": "FBA"}[analysis_id]
        whole = scope == "whole model"
        fn = (lambda: fba.run_pfba(model)) if analysis_id == "pfba" else (lambda: fba.run_fba(model))
        want_shadow = analysis_id == "shadow_prices"
        self._run_job(fn, lambda r: self._on_fba_done(method, r, whole, want_shadow, scope),
                      title=f"{method} · {scope}", kind=analysis_id)

    def _run_flux_sampling(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        from ..core.analysis import sampling
        from .dialogs.plot_gallery_dialog import PlotGalleryDialog
        from .viz import plots
        n = cfg.get("n_samples", 500)
        rxns = None
        if scope != "whole model" and self.project.categories.has(scope):
            rxns = list(self.project.categories.get(scope).reaction_ids)

        def done(samples):
            self._last_table = samples
            entries = [("Flux sampling — violin distributions",
                        lambda c: c.render(plots.flux_sampling_violin, samples,
                                           title="flux_sampling"))]
            self.status_label.setText(
                f"Flux sampling complete ({len(samples)} samples). Opening plot…")
            PlotGalleryDialog(self, entries).exec()
        self._run_job(lambda: sampling.run_flux_sampling(model, n_samples=n, reaction_ids=rxns),
                      done, title=f"Flux sampling ({n}) · {scope}", kind="flux_sampling")

    def _run_fva(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        def done(df):
            note = df.attrs.get("note", "")
            header = f"FVA · {scope}"
            if note:
                header += f" · ⚠ {note}"
            self._show_table(df.reset_index(), header)
        self._run_job(
            lambda: fba.run_fva(model, fraction_of_optimum=cfg["fraction_of_optimum"],
                                loopless=cfg["loopless"]),
            done, title=f"Flux Variability Analysis · {scope}", kind="fva")

    def _run_mutant(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        knockouts = self.explorer.selected_reaction_ids()
        if not knockouts:
            QMessageBox.information(
                self, "Mutant simulation",
                "Select one or more reactions in the Explorer (Reactions tab) to knock out, "
                "then run the mutant simulation.")
            return
        method = "room" if cfg["method"].startswith("ROOM") else "moma"
        whole = scope == "whole model"

        def done(result):
            self.analysis_panel.show_table(
                result.flux_table(),
                f"{result.method} mutant ({len(knockouts)} KO) · WT growth "
                f"{result.wt_growth:.4g} → mutant {result.growth:.4g} · "
                f"{result.n_changed} reactions changed")
            if whole and result.status == "optimal":
                fluxes = result.fluxes.to_dict()
                self.explorer.set_fluxes(fluxes)
                self.network_view.set_fluxes(fluxes)
            self.objective_label.setText(f"{result.method}: growth {result.growth:.4g}")
        self._run_job(lambda: mutant.run_mutant(model, knockouts, method=method), done,
                      title=f"Mutant simulation ({cfg['method'].split()[0]}) · {scope}",
                      kind="mutant")

    def _run_single_reaction_deletion(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        self._run_job(lambda: deletions.single_reaction_deletions(model),
                      lambda df: self._show_table(df, f"Single reaction deletion · {scope}"),
                      title=f"Single reaction deletion · {scope}",
                      kind="single_reaction_deletion")

    def _run_single_gene_deletion(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        self._run_job(lambda: deletions.single_gene_deletions(model),
                      lambda df: self._show_table(df, f"Single gene deletion · {scope}"),
                      title=f"Single gene deletion · {scope}", kind="single_gene_deletion")

    def _run_production_envelope(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        target = cfg["target"]
        if not model.reactions.has_id(target):
            QMessageBox.warning(self, "Production envelope", f"No reaction '{target}'.")
            return

        def done(df):
            self._show_table(df, f"Production envelope: {target} · {scope}")
            self.analysis_panel.plot_view.plot_production_envelope(df, target)
            self.analysis_panel.result_tabs.setCurrentWidget(self.analysis_panel.plot_view)
        self._run_job(lambda: phenotype.run_production_envelope(model, target, points=cfg["points"]), done,
                      title=f"Production envelope: {target} · {scope}", kind="production_envelope")

    def _run_robustness(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        control = cfg["control"]
        if not model.reactions.has_id(control):
            QMessageBox.warning(self, "Robustness", f"No reaction '{control}'.")
            return
        lower = None if cfg["auto_range"] else cfg["lower"]
        upper = None if cfg["auto_range"] else cfg["upper"]

        def done(df):
            self._show_table(df, f"Robustness: {control} · {scope}")
            self.analysis_panel.plot_view.plot_robustness(df, control)
            self.analysis_panel.result_tabs.setCurrentWidget(self.analysis_panel.plot_view)
        self._run_job(
            lambda: phenotype.run_robustness(model, control, points=cfg["points"],
                                             lower=lower, upper=upper), done,
            title=f"Robustness scan: {control} · {scope}", kind="robustness")

    def _run_phase_plane(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        rx, ry = cfg["reaction_x"], cfg["reaction_y"]
        if not (model.reactions.has_id(rx) and model.reactions.has_id(ry)):
            QMessageBox.warning(self, "Phase plane", "Both reactions must exist in the model.")
            return

        def done(df):
            self._show_table(df, f"Phase plane: {rx} vs {ry} · {scope}")
            self.analysis_panel.plot_view.plot_phase_plane(df, rx, ry)
            self.analysis_panel.result_tabs.setCurrentWidget(self.analysis_panel.plot_view)
        self._run_job(lambda: phenotype.run_phase_plane(model, rx, ry, points=cfg["points"]), done,
                      title=f"Phenotypic phase plane: {rx} vs {ry}", kind="phase_plane")

    # ----- Phase 6: community & pathways -------------------------------
    def build_community(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select two or more organism models", "", _MODEL_FILTER)
        if not paths:
            return
        if len(paths) < 2:
            QMessageBox.information(self, "Community model",
                                    "Select at least two model files to build a community.")
            return
        try:
            models = [io_models.load_model(p) for p in paths]
            names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
            cm = community_core.build_community(models, names)
        except (io_models.ModelLoadError, community_core.CommunityError) as exc:
            QMessageBox.critical(self, "Could not build community", str(exc))
            return
        self.community = cm
        self.project = Project.from_model(cm.model)
        self.project.settings["community_members"] = cm.member_biomass
        self._last_result = None
        self._refresh_all()
        self.status_label.setText(
            f"Built community of {len(cm.member_names)} members: {', '.join(cm.member_names)}.")
        QMessageBox.information(
            self, "Community model built",
            f"Combined {len(cm.member_names)} organisms into one community model "
            f"({len(cm.model.reactions)} reactions). Run 'Community member growth' in the "
            "Analysis tab to see each organism's growth. Use <b>Consortia objective "
            "(dominance &amp; min growth)</b> in the objective settings to control which "
            "species dominates and to stop the faster grower starving the others.")

    def _load_pathway_database(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open universal reaction database", "", _MODEL_FILTER)
        if not path:
            return
        try:
            db = io_models.load_model(path)
        except io_models.ModelLoadError as exc:
            QMessageBox.critical(self, "Could not load database", str(exc))
            return
        # a loaded local file is not part of our cache, so it is not deletable
        self._add_database(os.path.basename(path), db, cache_path=None, source="file")

    # ----- reaction-database registry ----------------------------------
    # A database entry: {name, source, cache_path, loader (callable|None),
    #   model (None until loaded), selected (checkbox), reactions, metabolites}.
    # Loading and the expensive combined-universal merge are deferred: listing is
    # instant, loading happens on "Load selected", merging happens at Predict.

    def _pathway_cache_paths(self) -> dict:
        enrich = databases.metanetx_reference_available()
        d = cache.databases_dir()
        return {
            "bigg": os.path.join(d, "bigg_universal_model.json"),
            "mnx": os.path.join(d, f"metanetx_universal_enzymatic{'_named2' if enrich else ''}.json"),
            "seed": os.path.join(d, "modelseed_universal.json"),
        }

    # ---- persistent registry of downloaded databases (fix 4) -----------
    def _registry_path(self) -> str:
        return os.path.join(cache.databases_dir(), "db_registry.json")

    def _load_registry(self) -> dict:
        import json
        try:
            with open(self._registry_path(), "r", encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}

    def _record_persistent_db(self, name, cache_path, source, reactions, metabolites) -> None:
        """Remember a downloaded/fetched database so it can be re-loaded offline
        on a later launch (any DB with a real cached file on disk)."""
        if not cache_path or not os.path.exists(cache_path):
            return
        import json
        reg = self._load_registry()
        reg[name] = {"cache_path": cache_path, "source": source or "fetch",
                     "reactions": int(reactions or 0), "metabolites": int(metabolites or 0)}
        try:
            with open(self._registry_path(), "w", encoding="utf-8") as fh:
                json.dump(reg, fh, indent=2)
        except OSError:
            pass

    def _registered_db_specs(self) -> list:
        """(name, source, cache_path, loader, reactions_hint) for every previously
        downloaded database still present on disk — available offline."""
        specs = []
        for name, info in self._load_registry().items():
            cpath = info.get("cache_path")
            if not cpath or not os.path.exists(cpath):
                continue

            def _loader(p=cpath):
                model = io_models.load_model(p)
                return pathway_design.make_model_ids_readable(model)
            specs.append((name, info.get("source", "offline cache"), cpath, _loader,
                          info.get("reactions", 0)))
        return specs

    #: Files in the databases folder that are caches, chemistry lookups or bundled example
    #: organism models — not reaction databases anyone would search against.
    _NON_DATABASE_FILES = {
        "db_registry.json",
        "bigg_metabolite_chem.json",
        "bigg_e_coli_core.json",
    }

    def _discovered_db_specs(self) -> list:
        """Model files sitting in the databases folder that nothing else lists.

        The registry is a convenience, not the source of truth. It is a single small JSON
        file: if it is deleted, never written, or lost when the folder is copied between
        machines, a merged database that cost hours to build silently disappears from the
        list while still occupying 25 MB on disk. Scanning the folder makes the list
        self-healing — the databases you have are the databases you are offered.
        """
        directory = cache.databases_dir()
        try:
            entries = sorted(os.listdir(directory))
        except OSError:
            return []

        known = {os.path.normcase(p) for p in self._pathway_cache_paths().values()}
        known |= {os.path.normcase(info.get("cache_path", ""))
                  for info in self._load_registry().values()}
        # Each default database is built through intermediates that share its stem —
        # `metanetx_universal_enzymatic`, `..._named`, `..._named2`. They are the same
        # data at different stages, so listing all three would offer the user a choice
        # between a database and two older copies of it.
        default_stems = tuple(
            os.path.normcase(os.path.splitext(os.path.basename(p))[0])
            for p in self._pathway_cache_paths().values() if p)

        specs = []
        for filename in entries:
            if not filename.endswith(".json") or filename in self._NON_DATABASE_FILES:
                continue
            path = os.path.join(directory, filename)
            stem = os.path.normcase(os.path.splitext(filename)[0])
            # Matched both ways: a build intermediate can be either a longer name than
            # the default in use (`..._named2`) or a shorter one (`..._enzymatic`),
            # depending on which stage the default currently points at.
            superseded = any(stem.startswith(d) or d.startswith(stem)
                             for d in default_stems)
            if os.path.normcase(path) in known or superseded:
                continue
            if not self._looks_like_model(path):
                continue

            def _loader(p=path):
                return pathway_design.make_model_ids_readable(io_models.load_model(p))

            specs.append((self._database_display_name(filename), "found on disk", path,
                          _loader, 0))
        return specs

    @staticmethod
    def _looks_like_model(path: str) -> bool:
        """Cheap check that a JSON file is a cobra model, without parsing megabytes.

        Deliberately no size floor. A single fetched pathway can be well under a
        kilobyte and is still a database the user asked for and expects to see.
        """
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                head = fh.read(65536)
        except OSError:
            return False
        return '"metabolites"' in head or '"reactions"' in head

    @staticmethod
    def _database_display_name(filename: str) -> str:
        stem = os.path.splitext(filename)[0].replace("_", " ").strip()
        return (stem[:1].upper() + stem[1:]) if stem else filename

    def _default_database_specs(self) -> list:
        """(name, source, cache_path, loader, reactions_hint) for the standard DBs."""
        p = self._pathway_cache_paths()

        def _bigg():
            return pathway_design.load_readable_cached(
                p["bigg"], pathway_design.download_bigg_universal)

        def _mnx():
            return pathway_design.load_readable_cached(
                p["mnx"], lambda: databases.build_metanetx_universal(
                    only_balanced=True, only_enzymatic=True,
                    enrich=databases.metanetx_reference_available()))

        def _seed():
            return pathway_design.load_readable_cached(
                p["seed"], databases.build_modelseed_universal)

        specs = []
        # Bundled offline universal first — always available, no network (Issue 4).
        if pathway_design.bundled_universal_available():
            specs.append(("Offline universal (bundled)", "offline", None,
                          pathway_design.load_bundled_universal, 15480))
        specs += [
            ("BiGG universal", "default", p["bigg"], _bigg, 28301),
            ("MetaNetX universal (enzymatic)", "default", p["mnx"], _mnx, 17495),
            ("ModelSEED universal", "default", p["seed"], _seed, 43000),
        ]
        return specs

    def _populate_available_databases(self) -> None:
        """List all databases the user can load — defaults plus any already loaded —
        unchecked and WITHOUT loading any model (instant)."""
        existing = {e["name"]: e for e in self._databases}
        dbs, seen = [], set()
        # Defaults, everything the registry remembers, and anything else on disk that
        # looks like a database — the last of those is what makes a merge survive a lost
        # registry.
        for name, source, cpath, loader, rxn_hint in (
                self._default_database_specs() + self._registered_db_specs()
                + self._discovered_db_specs()):
            if name in seen:
                continue
            seen.add(name)
            if name in existing:
                dbs.append(existing[name])
            else:
                dbs.append({"name": name, "source": source, "cache_path": cpath,
                            "loader": loader, "model": None, "selected": False,
                            "reactions": rxn_hint, "metabolites": 0})
        for name, e in existing.items():          # keep loaded local-file / fetched DBs
            if name not in seen:
                dbs.append(e)
        self._databases = dbs
        self._update_databases()

    def _add_database(self, name: str, model, *, cache_path=None, source: str = "") -> None:
        """Register an already-loaded/fetched database (model in hand), selected."""
        self._databases = [e for e in self._databases if e["name"] != name]
        self._databases.append({
            "name": name, "model": model, "selected": True, "loader": None,
            "cache_path": cache_path, "source": source,
            "reactions": len(model.reactions), "metabolites": len(model.metabolites),
        })
        self._pathway_db = None   # invalidate the cached combined universal
        # Persist to the registry so it's offline-loadable on a later launch (fix 4).
        self._record_persistent_db(name, cache_path, source,
                                   len(model.reactions), len(model.metabolites))
        self._update_databases(f"Loaded {name} ({len(model.reactions)} reactions).")

    def _load_selected_databases(self) -> None:
        """Reconcile loaded databases with the current selection: load ticked ones
        that aren't in memory yet, and *unload* unticked ones (freeing memory)."""
        # Unload databases that are no longer ticked but still hold a model.
        unloaded = 0
        for e in self._databases:
            if not e["selected"] and e["model"] is not None and e.get("loader") is not None:
                e["model"] = None
                unloaded += 1
        todo = [e for e in self._databases if e["selected"] and e["model"] is None
                and e["loader"] is not None]
        if not todo:
            self._pathway_db = None
            msg = f"Unloaded {unloaded} database(s)." if unloaded else None
            if not unloaded and not any(e["selected"] for e in self._databases):
                QMessageBox.information(
                    self, "Load databases",
                    "Tick one or more databases in the list first, then click "
                    "“Load selected databases”.")
            self._update_databases(msg or "")
            self._refresh_pathway_targets()
            return
        names = [e["name"] for e in todo]

        def work():
            out = []
            for e in todo:
                out.append((e["name"], e["loader"]()))
            return out

        def finish(loaded):
            for name, model in loaded:
                for e in self._databases:
                    if e["name"] == name and model is not None:
                        e["model"] = model
                        e["reactions"] = len(model.reactions)
                        e["metabolites"] = len(model.metabolites)
            self._pathway_db = None
            tail = f" · unloaded {unloaded}" if unloaded else ""
            self._update_databases(f"Loaded {len(loaded)} database(s){tail}.")
            self._refresh_pathway_targets()

        ok, res = run_busy(
            self, f"Loading {len(names)} database(s): {', '.join(names)}…",
            work, title="Loading reaction databases",
            after=finish, after_message="Preparing databases…", cancelable=True)
        if not ok and not was_cancelled(res):
            self._show_db_load_error(res)

    def _rename_database(self, name: str) -> None:
        """Give a database a friendlier display name.

        Merges get an accurate but unwieldy auto-name ("BiGG + MetaNetX + ModelSEED
        (merged)"); the user should be able to call it whatever they like.
        """
        entry = next((e for e in self._databases if e["name"] == name), None)
        if entry is None:
            return
        new, ok = QInputDialog.getText(self, "Rename database", "Display name:",
                                       text=name)
        new = (new or "").strip()
        if not ok or not new or new == name:
            return
        if any(e["name"] == new for e in self._databases):
            QMessageBox.information(self, "Rename database",
                                    f"A database called “{new}” already exists — "
                                    f"pick another name.")
            return
        entry["name"] = new
        # The name is the key used for selection/removal elsewhere, so anything that
        # remembers it by name has to move with it.
        if getattr(self, "_pathway_db_names", None) and name in self._pathway_db_names:
            self._pathway_db_names = [new if n == name else n
                                      for n in self._pathway_db_names]
        self._update_databases(f"Renamed database to “{new}”.")

    def _merge_databases(self) -> None:
        """Unify all loaded databases into one deduplicated cytosolic database (#B6):
        each compound/reaction appears exactly once (no more isoprene ×4)."""
        included = [(e["name"], e["model"]) for e in self._databases
                    if e["model"] is not None]
        if len(included) < 1:
            QMessageBox.information(self, "Merge databases",
                                   "Load one or more databases first (tick them and click "
                                   "“Load selected databases”).")
            return

        # Merging unifies every compound across databases by cross-reference, which is
        # inherently expensive: merging BiGG + MetaNetX + SEED (~71k reactions) took
        # nearly three hours in testing. Say so BEFORE the user commits to it.
        total_rxn = sum(len(m.reactions) for _n, m in included)
        est = "a few minutes"
        if total_rxn > 50000:
            est = "up to several HOURS"
        elif total_rxn > 20000:
            est = "many minutes, possibly over an hour"
        names = ", ".join(n for n, _m in included)
        box = QMessageBox(self)
        box.setWindowTitle("Merge databases — this takes a long time")
        box.setIcon(QMessageBox.Warning)
        box.setTextFormat(Qt.RichText)
        box.setText(
            f"<b>Merging is a long, one-off computation.</b><br><br>"
            f"About to merge <b>{len(included)}</b> database(s) — {names} — "
            f"totalling <b>{total_rxn:,}</b> reactions.<br><br>"
            f"Estimated time: <b>{est}</b>. Every compound is matched against every "
            f"other by cross-reference, so the cost grows quickly with size.<br><br>"
            f"The result is <b>saved as a new database</b>, so this only has to be done "
            f"once — afterwards you can load it instantly from the database list. "
            f"Your existing databases are kept.<br><br>Proceed?")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.Cancel)
        box.setDefaultButton(QMessageBox.Cancel)
        if box.exec() != QMessageBox.Yes:
            return

        from ..core import pathway_search

        merged_name = " + ".join(n for n, _m in included) + " (merged)"

        def work():
            merged = pathway_search.merge_databases([m for _n, m in included])
            # Persist immediately: this cost hours, so it must survive the session and
            # be loadable like any other database.
            path = None
            try:
                import re

                from ..core.cache import databases_dir
                slug = re.sub(r"[^A-Za-z0-9]+", "_", merged_name).strip("_")[:60]
                path = os.path.join(databases_dir(), f"{slug}.json")
                io_models.save_model(merged, path)
            except Exception:  # noqa: BLE001 — a save failure must not lose the merge
                path = None
            return merged, path

        def finish(built):
            merged, path = built
            # Keep the source databases: the merge is an ADDITIONAL database, not a
            # replacement. Destroying the originals would make the merge irreversible
            # and force a several-hour rebuild to get back. Untick them so the search
            # uses the merge alone (combining merge + sources would double everything).
            for e in self._databases:
                e["selected"] = False
            # _add_database registers it in the on-disk registry too, which is what
            # makes a merge that cost hours show up in the list on the next launch.
            self._add_database(merged_name, merged, cache_path=path, source="merged")
            self._refresh_pathway_targets()
            if path:
                QMessageBox.information(
                    self, "Merge complete",
                    f"Merged {len(included)} database(s) into "
                    f"“{merged_name}”:\n\n"
                    f"  {len(merged.reactions):,} reactions\n"
                    f"  {len(merged.metabolites):,} unique metabolites\n\n"
                    f"It has been saved and added to your database list, so you will "
                    f"not have to merge again — it loads instantly next time.\n\n"
                    f"You can rename it from “Manage reaction databases…”.")
            else:
                QMessageBox.warning(
                    self, "Merge complete (not saved)",
                    f"The merge succeeded ({len(merged.reactions):,} reactions) but "
                    f"could not be written to disk, so it exists only in this session. "
                    f"Use “Save loaded database…” to keep it.")

        ok, res = run_busy(self, f"Merging {total_rxn:,} reactions (removing duplicates) — "
                                 f"this can take {est}…", work,
                           title="Merge databases", after=finish,
                           after_message="Rebuilding target list…", cancelable=True)
        if not ok and not was_cancelled(res):
            QMessageBox.critical(self, "Could not merge databases", str(res))

    def _save_loaded_database(self) -> None:
        """Save the currently loaded/merged database to a JSON file for offline reuse
        (loadable later via “Load reaction database…”) (#P2)."""
        combined = self._combined_included_model()
        if combined is None:
            QMessageBox.information(self, "Save loaded database",
                                   "No database is loaded. Tick and load databases first "
                                   "(and optionally merge them).")
            return
        from .widgets.dialog_util import choose_save_path
        path = choose_save_path(self, "Save loaded database", "reaction_database.json",
                                "Reaction database (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        ok, res = run_busy(
            self, f"Saving database ({len(combined.reactions)} reactions)…",
            lambda: io_models.save_model(combined, path), title="Save database",
            cancelable=True)
        if not ok:
            if not was_cancelled(res):
                QMessageBox.critical(self, "Could not save database", str(res))
            return
        self.status_label.setText(f"Saved database to {os.path.basename(path)}.")

    def _combined_included_model(self):
        """A single model combining all selected+loaded databases (union of reactions).

        Cached in ``self._pathway_db``; rebuilt only when the selection changes. This
        merge is expensive on genome-scale universals, so it is built lazily (at
        prediction time), never during listing/loading."""
        if self._pathway_db is not None:
            return self._pathway_db
        included = [e["model"] for e in self._databases
                    if e["selected"] and e["model"] is not None]
        if not included:
            return None
        if len(included) == 1:
            self._pathway_db = included[0]
            return self._pathway_db
        merged = cobra.Model("combined_databases")
        seen = set()
        for m in included:
            new = []
            for rxn in m.reactions:
                if rxn.id in seen:
                    continue
                seen.add(rxn.id)
                new.append(rxn.copy())
            merged.add_reactions(new)
        self._pathway_db = merged
        return merged

    def _update_databases(self, status: str = "") -> None:
        """Refresh the panel summary only — cheap; no merge, no target rebuild."""
        entries = [{"name": e["name"], "reactions": e["reactions"],
                    "metabolites": e.get("metabolites", 0),
                    "selected": e["selected"], "loaded": e["model"] is not None,
                    "source": e.get("source", "")}
                   for e in self._databases]
        self.pathway_panel.set_databases_summary(entries)
        if status:
            self.status_label.setText(status)

    def _set_database_included(self, name: str, selected: bool) -> None:
        for e in self._databases:
            if e["name"] == name:
                e["selected"] = selected
        self._pathway_db = None   # selection changed -> rebuild combined next predict
        self._update_databases()

    def _remove_database(self, name: str) -> None:
        entry = next((e for e in self._databases if e["name"] == name), None)
        if entry is None:
            return
        path = entry.get("cache_path")
        if path and os.path.exists(path):
            try:
                os.remove(path)
            except OSError:
                pass
        # Drop the cache pickle too so it is not silently reused.
        pkl = (path + ".readable.pkl") if path else None
        if pkl and os.path.exists(pkl):
            try:
                os.remove(pkl)
            except OSError:
                pass
        self._databases = [e for e in self._databases if e["name"] != name]
        self._pathway_db = None
        self._update_databases(f"Removed database {name}.")

    _DB_SOURCE_LABELS = {"fetch": "downloaded", "file": "local file", "demo": "demo"}

    def _reconcile_identifiers(self) -> None:
        """Propose, review and store metabolite correspondences with the loaded database.

        Most compounds already resolve through their cross-references, so this covers the
        residue: entries whose annotation sets do not overlap at all. Left unresolved they
        enter the model as orphans and take a route's flux to zero without any visible
        error, which is the failure this exists to prevent.
        """
        if self.project is None:
            return
        combined = self._combined_included_model()
        if combined is None:
            QMessageBox.information(
                self, "Reconcile metabolite identifiers",
                "Load a reaction database first, in the Pathway Design tab.")
            return

        from ..core import id_reconcile as rec
        from .dialogs.reconcile_dialog import ReconcileDialog

        stored = self.project.datasets.get(rec.PROJECT_KEY)
        preselected = rec.from_record(stored)

        def work():
            return rec.propose(self.project.model, combined)

        ok, proposals = run_busy(
            self, "Comparing identifiers between the model and the database…", work,
            title="Reconcile metabolite identifiers", cancelable=True)
        if not ok:
            if not was_cancelled(proposals):
                QMessageBox.critical(self, "Could not compare identifiers",
                                     str(proposals))
            return
        if not proposals:
            QMessageBox.information(
                self, "Reconcile metabolite identifiers",
                "Nothing left to reconcile. Every database metabolite this model needs "
                "already resolves through its cross-references.")
            return

        dialog = ReconcileDialog(self, proposals,
                                 database_name=", ".join(
                                     e["name"] for e in self._databases
                                     if e.get("selected") and e.get("model") is not None),
                                 preselected=preselected)
        if dialog.exec() != QDialog.Accepted:
            return

        record = dialog.record()
        self.project.datasets[rec.PROJECT_KEY] = record
        self.project._dirty = True
        count = len(record.get("entries", []))
        QMessageBox.information(
            self, "Mapping saved",
            f"{count} correspondence(s) stored with the project. They are listed in any "
            f"report this project produces, and can be changed from this dialog at any "
            f"time.")

    def _manage_databases(self) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle("Manage reaction databases")
        # Fit neatly within the interface: cap the width to a fraction of the screen.
        screen = self.screen().availableGeometry() if self.screen() else None
        max_w = min(600, int(screen.width() * 0.6)) if screen else 600
        dlg.setMinimumWidth(420)
        dlg.setMaximumWidth(max_w)
        dlg.resize(max_w, 420)
        if not self._databases:
            self._populate_available_databases()
        outer = QVBoxLayout(dlg)
        outer.addWidget(QLabel(
            "Reaction databases available for pathway design. Tick the ones you want and click "
            "“Load selected databases” to load them into memory (ticking alone doesn’t load). "
            "Loaded, ticked databases are combined for the search. The trash button permanently "
            "deletes a downloaded database from disk; local-file databases are only unloaded."))

        host = QWidget()
        vbox = QVBoxLayout(host)
        vbox.setAlignment(Qt.AlignTop)

        def rebuild():
            while vbox.count():
                it = vbox.takeAt(0)
                w = it.widget()
                if w is not None:
                    w.deleteLater()
            if not self._databases:
                empty = QLabel("No databases available.")
                empty.setWordWrap(True)
                vbox.addWidget(empty)
            for e in self._databases:
                row = QWidget()
                h = QHBoxLayout(row)
                h.setContentsMargins(2, 2, 2, 2)
                src = self._DB_SOURCE_LABELS.get(e.get("source", ""), e.get("source", ""))
                mets = e.get("metabolites", 0)
                loaded = e.get("model") is not None
                col = QVBoxLayout()
                col.setSpacing(0)
                chk = QCheckBox(e["name"])
                chk.setChecked(e["selected"])
                chk.toggled.connect(lambda on, n=e["name"]: self._set_database_included(n, on))
                state = "loaded" if loaded else "not loaded"
                count = (f"{e['reactions']:,} reactions · {mets:,} metabolites" if loaded
                         else f"~{e['reactions']:,} reactions")
                detail = QLabel(f"{count} · {state}" + (f" · {src}" if src else ""))
                detail.setStyleSheet("color:#5f6368; font-size:11px; padding-left:20px;")
                col.addWidget(chk)
                col.addWidget(detail)
                rename = QToolButton()
                rename.setText("Rename")
                rename.setToolTip("Give this database a name that means something to you "
                                  "(e.g. shorten a long auto-generated merge name).")
                rename.setAutoRaise(True)
                rename.clicked.connect(
                    lambda _=False, n=e["name"]: (self._rename_database(n), rebuild()))
                trash = QToolButton()
                trash.setText("Remove")
                trash.setToolTip("Delete this database from permanent storage")
                trash.setAutoRaise(True)
                trash.clicked.connect(
                    lambda _=False, n=e["name"]: (self._confirm_remove_db(n), rebuild()))
                h.addLayout(col, 1)
                h.addWidget(rename)
                h.addWidget(trash)
                vbox.addWidget(row)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(host)
        outer.addWidget(scroll, 1)

        # Bulk actions
        actions = QHBoxLayout()
        inc_all = QPushButton("Select all")
        inc_all.clicked.connect(lambda: (self._set_all_included(True), rebuild()))
        exc_all = QPushButton("Deselect all")
        exc_all.clicked.connect(lambda: (self._set_all_included(False), rebuild()))
        load_sel = QPushButton("Load selected databases")
        load_sel.setObjectName("primary")
        load_sel.setToolTip("Download (if needed) and load the ticked databases into memory.")
        load_sel.clicked.connect(lambda: (dlg.accept(), self._load_selected_databases()))
        actions.addWidget(inc_all)
        actions.addWidget(exc_all)
        actions.addStretch(1)
        actions.addWidget(load_sel)
        outer.addLayout(actions)

        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        outer.addWidget(buttons)
        rebuild()
        dlg.exec()

    def _set_all_included(self, included: bool) -> None:
        for e in self._databases:
            e["selected"] = included
        self._pathway_db = None
        self._update_databases()

    def _confirm_remove_db(self, name: str) -> None:
        entry = next((e for e in self._databases if e["name"] == name), None)
        if entry is None:
            return
        deletable = bool(entry.get("cache_path"))
        msg = (f"Permanently delete “{name}” from disk?" if deletable else
               f"Unload “{name}”? (It was loaded from a local file and won't be deleted.)")
        if QMessageBox.question(self, "Remove database", msg) == QMessageBox.Yes:
            self._remove_database(name)

    # BiGG universal + MetaNetX are loaded automatically on startup, so they are
    # not offered here; these are the *additional* on-demand sources.
    _PATHWAY_SOURCES = [
        "A specific BiGG model by id (e.g. iML1515)…",
        "ModelSEED universal biochemistry (SEED namespace, ~40k reactions)",
        "KEGG — focused database around a target product…",
        "MetaNetX universal (rebuild / with names — advanced)",
    ]

    def _fetch_pathway_database(self) -> None:
        choice, ok = QInputDialog.getItem(
            self, "Fetch reaction database online",
            "Source:", self._PATHWAY_SOURCES, 0, False)
        if not ok:
            return
        idx = self._PATHWAY_SOURCES.index(choice)

        if idx == 0:
            model_id, ok = QInputDialog.getText(self, "BiGG model id", "BiGG model id:")
            if ok and model_id.strip():
                self._fetch_bigg(universal=False, model_id=model_id.strip())
        elif idx == 1:
            self._fetch_modelseed()
        elif idx == 2:
            self._fetch_kegg()
        elif idx == 3:
            self._fetch_metanetx()

    def _run_fetch(self, message: str, work, *, then=None) -> None:
        """Run a database download/build on a background thread behind a busy
        dialog and register the result. Downloads stay *in-process* (they are not
        shipped back through a job queue) so large genome-scale universals load
        reliably. ``work`` returns ``(model, name, cache_path)`` or a list of such
        tuples."""
        ok, result = run_busy(self, message, work, title="Loading reaction database",
                              cancelable=True)
        if not ok:
            if not was_cancelled(result):
                self._show_db_load_error(result)
            return
        items = result if isinstance(result, list) else [result]
        try:
            for model, name, cache_path in items:
                self._add_database(name, model, cache_path=cache_path, source="fetch")
        except Exception as exc:  # noqa: BLE001 - surface instead of silently dropping
            QMessageBox.critical(self, "Could not register database", str(exc))
            return
        if then:
            then()

    def _show_db_load_error(self, err) -> None:
        """Clear, actionable message on a failed database download (Issue 4)."""
        text = str(err)
        looks_offline = any(s in text for s in ("403", "timed out", "timeout", "URLError",
                                                "Failed to establish", "getaddrinfo",
                                                "Download failed", "Name or service"))
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("Could not load database")
        if looks_offline:
            box.setText("The reaction database could not be downloaded — you may be "
                        "offline or behind a proxy.")
            hint = ""
            if pathway_design.bundled_universal_available():
                hint = ("• Tick “Offline universal (bundled)” in the Reaction databases "
                        "panel and click “Load selected databases” — it needs no internet.\n")
            box.setInformativeText(
                hint +
                "• Or use “Load database from file…” to load a universal model you "
                "already have (SBML/JSON).\n"
                "• Or retry when a connection is available.")
        else:
            box.setText("The reaction database could not be loaded.")
            box.setInformativeText(text)
        box.exec()

    def _load_default_databases(self, *, then=None, prompt: bool = True) -> None:
        """Select and load the recommended broad databases (BiGG universal +
        MetaNetX enzymatic). Just marks them selected and routes through the normal
        "load selected" path (which downloads once, caches, and loads fast after)."""
        if not self._databases:
            self._populate_available_databases()
        defaults = {"BiGG universal", "MetaNetX universal (enzymatic)"}
        for e in self._databases:
            if e["name"] in defaults:
                e["selected"] = True
        self._pathway_db = None
        self._update_databases()
        self._load_selected_databases()
        if then:
            then()

    def _fetch_bigg(self, *, universal: bool, model_id: str = "") -> None:
        if universal:
            confirm = QMessageBox.question(
                self, "Download database",
                "Download the BiGG universal model (tens of MB) from bigg.ucsd.edu? "
                "It is cached for reuse.")
            if confirm != QMessageBox.Yes:
                return
        fname = "bigg_universal_model.json" if universal else f"bigg_{model_id}.json"
        cache_path = os.path.join(cache.databases_dir(), fname)

        def work():
            m = (pathway_design.download_bigg_universal() if universal
                 else pathway_design.fetch_bigg_model(model_id))
            pathway_design.make_model_ids_readable(m)
            return m, ("BiGG universal" if universal else f"BiGG {model_id}"), cache_path
        self._run_fetch(
            f"Downloading BiGG {'universal model' if universal else model_id}…", work)

    def _fetch_modelseed(self) -> None:
        confirm = QMessageBox.question(
            self, "ModelSEED biochemistry",
            "Download and build the ModelSEED universal biochemistry database "
            "(reactions + compounds, ~15 MB from the ModelSEEDDatabase project, cached "
            "for reuse)?\n\nIt uses the SEED namespace and carries compound names and EC "
            "numbers, giving broad coverage complementary to BiGG and MetaNetX.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if confirm != QMessageBox.Yes:
            return
        cache_path = os.path.join(cache.databases_dir(), "modelseed_universal.json")

        def work():
            db = databases.build_modelseed_universal()
            pathway_design.make_model_ids_readable(db)
            return db, "ModelSEED universal", cache_path
        self._run_fetch("Downloading ModelSEED biochemistry…", work)

    def _fetch_metanetx(self) -> None:
        only_ec = QMessageBox.question(
            self, "MetaNetX universal database",
            "Build the MetaNetX universal reaction database (downloads a ~10 MB table "
            "from metanetx.org, cached for reuse).\n\n"
            "Restrict to enzyme-catalysed reactions only (those with an EC number)? "
            "This is smaller and focuses on realistic enzymatic steps; choose No to "
            "include all balanced reactions.") == QMessageBox.Yes

        # Compound names + cross-references come from a large one-time download.
        enrich = databases.metanetx_reference_available()
        if not enrich:
            enrich = QMessageBox.question(
                self, "MetaNetX names & cross-references",
                "Also download the MetaNetX reference data (compounds ~680 MB + reactions "
                "~80 MB, one-time, stored permanently)?\n\n"
                "With it, entries show their real names and readable ids, and can be "
                "matched to your model across identifier systems. Without it they "
                "appear as MNXM/MNXR ids and only match models that already carry "
                "MetaNetX or KEGG annotations.",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes) == QMessageBox.Yes

        downloading_ref = enrich and not databases.metanetx_reference_available()
        kind_tag = "enzymatic" if only_ec else "balanced"
        cache_path = os.path.join(
            cache.databases_dir(),
            f"metanetx_universal_{kind_tag}{'_named2' if enrich else ''}.json")

        def work():
            db = databases.build_metanetx_universal(
                only_balanced=True, only_enzymatic=only_ec, enrich=enrich)
            pathway_design.make_model_ids_readable(db)
            label = "MetaNetX universal" + (" (enzymatic)" if only_ec else "")
            return db, label, cache_path
        msg = ("Downloading MetaNetX reference data (~680 MB)…" if downloading_ref
               else "Building MetaNetX universal database…")
        self._run_fetch(msg, work)

    def _fetch_kegg(self) -> None:
        target, ok = QInputDialog.getText(
            self, "KEGG — focused database",
            "Target product (compound name or KEGG compound id):")
        if not ok or not target.strip():
            return
        steps, ok = QInputDialog.getInt(
            self, "KEGG — neighbourhood", "Expansion steps around the target "
            "(0 = only reactions touching the target; 1-2 = include neighbours):",
            1, 0, 3)
        if not ok:
            return

        def work():
            model, label, cache_path = databases.build_kegg_pathway_db(
                target.strip(), expand_steps=steps)
            return model, label, cache_path
        self._run_fetch(f"Fetching KEGG reactions around “{target.strip()}”…", work)

    @staticmethod
    def _target_display(name: str, met_id: str) -> str:
        base = met_id.rsplit("_", 1)[0]
        if name and name not in (met_id, base):
            return f"{name} ({met_id})"
        return met_id

    def _refresh_pathway_targets(self) -> None:
        """Build the target-metabolite list: host metabolites plus the metabolites of
        the selected, *loaded* databases (deduped cheaply by id/name). This iterates
        the already-loaded DB models directly — no merge — and is invoked behind the
        load dialog, so it never freezes the UI after loading finishes."""
        if self.project is None:
            self.pathway_panel.set_targets([])
            return
        host = self.project.model
        items = [(self._target_display(m.name, m.id), m.id) for m in host.metabolites]
        seen = {m.id for m in host.metabolites}
        host_names = {(m.name or "").strip().lower() for m in host.metabolites}
        host_names.discard("")
        for e in self._databases:
            if not (e.get("selected") and e.get("model") is not None):
                continue
            for m in e["model"].metabolites:
                if m.id in seen:
                    continue
                nm = (m.name or "").strip().lower()
                if nm and nm in host_names:
                    continue
                seen.add(m.id)
                items.append((self._target_display(m.name, m.id), m.id))
        items.sort(key=lambda it: it[0].lower())
        self.pathway_panel.set_targets(items)

    def _lookup_enzymes(self, ec_numbers: list, reaction_id: str) -> None:
        if not ec_numbers:
            return
        ec = ec_numbers[0]
        self.status_label.setText(f"Looking up UniProt enzymes for EC {ec}…")

        def fetch():
            return ec, databases.uniprot_enzymes_for_ec(ec, reviewed_only=True, limit=50)

        def done(result):
            ec_used, df = result
            self.status_label.setText(f"UniProt: {len(df)} enzyme(s) for EC {ec_used}.")
            self._show_enzyme_dialog(ec_used, reaction_id, df)
        self._run_job(fetch, done, title=f"UniProt enzyme lookup (EC {ec})…", kind="uniprot")

    def build_enzyme_dialog(self, ec: str, reaction_id: str, df) -> QDialog:
        """Build (but do not show) the UniProt enzyme-results dialog."""
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Enzymes for EC {ec}")
        dlg.resize(720, 440)
        lay = QVBoxLayout(dlg)
        header = QLabel(
            f"<b>Reviewed UniProt enzymes for EC {ec}</b> "
            f"(reaction {reaction_id}).<br>Candidate genes to express for this step.")
        header.setWordWrap(True)
        lay.addWidget(header)
        view = ResultsView()
        if df is None or df.empty:
            header.setText(header.text() + "<br><br>No reviewed UniProt entries found "
                           "for this EC number.")
        else:
            view.show_dataframe(df, f"UniProt enzymes — EC {ec}")
        lay.addWidget(view, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        return dlg

    def _show_enzyme_dialog(self, ec: str, reaction_id: str, df) -> None:
        self.build_enzyme_dialog(ec, reaction_id, df).exec()

    def _has_loaded_db(self) -> bool:
        return any(e.get("selected") and e.get("model") is not None for e in self._databases)

    def _resolve_target_in(self, combined, target: str) -> str:
        """Resolve a typed/selected target to a metabolite id in host or combined DB,
        falling back to a case-insensitive name match against the combined database."""
        host = self.project.model
        if host.metabolites.has_id(target) or (
                combined is not None and combined.metabolites.has_id(target)):
            return target
        low = target.strip().lower()
        if combined is not None and low:
            for m in combined.metabolites:            # exact name first
                if (m.name or "").strip().lower() == low:
                    return m.id
            for m in combined.metabolites:            # then contains
                if low in (m.name or "").strip().lower():
                    return m.id
        return target

    @staticmethod
    def _fetch_smiles_online(met) -> str:
        """Best-effort SMILES for a metabolite from online sources, via its
        cross-references: PubChem (by InChIKey, then name) and KEGG (molblock → SMILES).
        Returns "" if nothing resolves."""
        from .widgets import structure_fetcher as sf

        smiles, inchi, inchikey, kegg, _chebi = sf.metabolite_structure_hints(met)
        if smiles:
            return smiles
        name = getattr(met, "name", "") or ""
        try:
            s = sf._pubchem_smiles(inchikey=inchikey,
                                   name="" if sf._ambiguous_name(name) else name)
            if s:
                return s
        except Exception:  # noqa: BLE001
            pass
        try:
            mol = sf._kegg_molblock(kegg) if kegg else ""
            if mol:
                from rdkit import Chem
                m = Chem.MolFromMolBlock(mol, sanitize=True, removeHs=True)
                if m is not None:
                    return Chem.MolToSmiles(m)
        except Exception:  # noqa: BLE001
            pass
        if inchi:
            try:
                from rdkit import Chem
                m = Chem.MolFromInchi(inchi)
                if m is not None:
                    return Chem.MolToSmiles(m)
            except Exception:  # noqa: BLE001
                pass
        return ""

    def _retrorules_suggest(self, target: str) -> None:
        """Rule-based retrosynthesis (RetroRules) for the target metabolite.

        Unlike the database search, this proposes steps from generalised reaction RULES,
        so it can suggest chemistry no loaded database contains. It needs a structure
        (SMILES/InChI) for the target and the host's native InChIKeys; the result is a
        set of PREDICTED steps to verify, shown distinctly from database reactions.
        """
        if self.project is None:
            return
        from ..core import retrorules

        if not retrorules.rdkit_available():
            QMessageBox.information(
                self, "RetroRules unavailable",
                "Rule-based retrosynthesis needs RDKit, which is not available in this "
                "build.")
            return

        # Resolve the target to a structure. Prefer a database metabolite's SMILES/InChI
        # annotation; if absent (most genome-scale metabolites carry only an InChIKey or
        # a KEGG id), fetch the structure ONLINE (PubChem by InChIKey/name, then KEGG);
        # only ask the user as a last resort.
        combined = self._combined_included_model()
        resolved = self._resolve_target_in(combined, target)
        target_met = None
        target_name = target
        for mdl in (self.project.model, combined):
            if mdl is not None and mdl.metabolites.has_id(resolved):
                target_met = mdl.metabolites.get_by_id(resolved)
                target_name = target_met.name or resolved
                break

        smiles = retrorules.metabolite_smiles(target_met) if target_met is not None else ""
        if not smiles and target_met is not None:
            ok, res = run_busy(
                self, f"Fetching the structure of {target_name} online…",
                lambda m=target_met: self._fetch_smiles_online(m),
                title="RetroRules", cancelable=True)
            if ok and res:
                smiles = res
        if not smiles:
            smiles, ok = QInputDialog.getText(
                self, "RetroRules — target structure",
                f"No SMILES/InChI could be found or fetched for “{target_name}”.\n"
                "Enter a SMILES for the target to run rule-based retrosynthesis:")
            if not ok or not smiles.strip():
                return
            smiles = smiles.strip()

        # Offer to install the full ruleset if it is not present (bundled curated rules
        # work, but with limited coverage).
        ruleset_note = "Using the bundled curated ruleset (limited coverage)."
        if retrorules.assets_present():
            ruleset_note = ("Using the full RetroRules RR02 dataset "
                            f"(diameter {retrorules.DEFAULT_DIAMETER}).")
        else:
            ans = QMessageBox.question(
                self, "RetroRules dataset",
                "The full RetroRules dataset (~43 MB download) is not installed.\n\n"
                "Download it now for full coverage? (Choose No to use the small "
                "bundled ruleset.)",
                QMessageBox.Yes | QMessageBox.No)
            if ans == QMessageBox.Yes:
                ok, res = run_busy(
                    self, "Downloading & extracting RetroRules (~43 MB)…",
                    lambda: retrorules.install_full_ruleset(), title="RetroRules",
                    cancelable=False)
                if ok and retrorules.assets_present():
                    ruleset_note = ("Using the full RetroRules RR02 dataset "
                                    f"(diameter {retrorules.DEFAULT_DIAMETER}).")
                elif not ok:
                    QMessageBox.warning(self, "RetroRules download failed",
                                        f"{res}\n\nFalling back to the bundled ruleset.")

        host = self.project.model

        # Expanded input options (#3.2): how many alternative routes, how deep, and how
        # to rank them.
        from .dialogs.retrorules_dialog import RetroRulesDialog, RetroRulesOptionsDialog
        opt = RetroRulesOptionsDialog(self, target_name)
        if not opt.exec():
            return
        cfg = opt.values()

        # Capture the seed the search will use, so the results dialog can quote the seed
        # that actually produced these routes rather than whatever is saved later.
        from ..core import preferences as _prefs
        rr_seed = _prefs.retrorules_seed()

        def work():
            return retrorules.suggest_routes_multi(
                smiles, host, n_alternatives=cfg["n_alternatives"],
                rank_by=cfg["rank_by"], max_steps=cfg["max_steps"], time_budget=120.0,
                seed=rr_seed)

        ok, routes = run_busy(
            self, f"Searching reaction rules for {target_name}… (this can take a "
                  f"minute)", work, title="RetroRules", cancelable=True)
        if not ok:
            if not was_cancelled(routes):
                QMessageBox.warning(self, "RetroRules failed", str(routes))
            return
        routes = routes or []

        # Resolve every route compound to a human-readable name: first from the loaded
        # databases (by InChIKey), then ONLINE (PubChem) for the rule-generated
        # intermediates the databases don't contain. Done once here, behind a progress
        # dialog, so the names are ready for both the display and the COBRA reactions.
        all_smiles = {smiles}
        for route in routes:
            for s in route.steps:
                all_smiles.add(s["product"])
                all_smiles.update(s["precursors"])
            all_smiles.update(route.terminal_precursors)

        local = self._smiles_name_index(combined)

        def resolve_names():
            names = {}
            for smi in all_smiles:
                nm = local(smi)
                if not nm:
                    from .widgets.structure_fetcher import name_from_smiles
                    nm = name_from_smiles(smi)      # online, cached
                if nm:
                    names[smi] = nm
            return names

        ok, names = run_busy(self, "Looking up compound names…", resolve_names,
                             title="RetroRules", cancelable=True)
        if not ok:
            names = {smi: local(smi) for smi in all_smiles if local(smi)}
        name_for = names.get

        rank_label = retrorules.RANK_LABELS.get(cfg["rank_by"], "")
        dlg = RetroRulesDialog(self, target_name, routes, target_smiles=smiles,
                               rank_label=rank_label, ruleset_note=ruleset_note,
                               name_for=name_for, seed=rr_seed)
        dlg.add_reactions_requested.connect(
            lambda d=dlg, s=smiles, n=target_name, nm=names:
            self._add_retrorules_pathway(d.selected_route(), s, n, nm,
                                         step_indices=d.selected_step_indices()))
        dlg.exec()

    def _smiles_name_index(self, combined):
        """A fast SMILES → name resolver from an InChIKey index over the host + loaded
        databases (no network). Empty string when the compound is not recognised."""
        index: dict = {}
        for mdl in (self.project.model, combined):
            if mdl is None:
                continue
            for m in mdl.metabolites:
                ann = getattr(m, "annotation", None) or {}
                v = ann.get("inchi_key") or ann.get("inchikey")
                if not v:
                    continue
                for k in (v if isinstance(v, (list, tuple)) else [v]):
                    blk = str(k).split("-", 1)[0]
                    if blk and blk not in index and (m.name and m.name != m.id):
                        index[blk] = m.name

        def name_for(smiles: str) -> str:
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return ""
                blk = Chem.MolToInchiKey(mol).split("-", 1)[0]
                return index.get(blk, "")
            except Exception:  # noqa: BLE001
                return ""
        return name_for

    def _add_retrorules_pathway(self, route, target_smiles, target_name, names,
                                *, step_indices=None) -> None:
        """Load a rule-based route into Pathway Design as a suggested pathway, so the
        normal Apply flow can add it to the model. ``names`` maps SMILES → human name
        (already resolved locally + online), so the COBRA metabolites carry real names.
        ``step_indices`` (if given) restricts the pathway to the steps the user ticked."""
        if self.project is None or route is None or not route.steps:
            return
        from ..core import pathway_design, retrorules

        # Restrict to the user-selected steps (#3.2 checkboxes). Terminal precursors are
        # recomputed as the selected steps' precursors that no selected step produces.
        if step_indices is not None and len(step_indices) != len(route.steps):
            chosen = [route.steps[i] for i in step_indices if 0 <= i < len(route.steps)]
            if not chosen:
                return
            produced = {s["product"] for s in chosen}
            terminals = [p for s in chosen for p in s["precursors"] if p not in produced]
            route = retrorules.RuleRoute(
                target=route.target, steps=chosen,
                terminal_precursors=list(dict.fromkeys(terminals)),
                complete=route.complete)

        names = dict(names or {})
        model, ids, target_mid = retrorules.build_suggested_model(
            route, target_smiles, target_name=target_name, names=names)
        # Register the suggested reactions as a database so Apply/Display/Draw all work
        # unchanged (they operate on the combined included model).
        self._add_database("RetroRules suggestions", model, source="retrorules")

        result = pathway_design._result_from_ids(
            self.project.model, model, target_mid or target_name, ids, rank=0)
        # A rule route's flux is measured on a synthetic demand for a novel compound, so
        # it saturates against a generic host bound and is NOT comparable with a database
        # route's flux. Mark it so the UI reports it qualitatively (L2).
        result.flux_is_indicative = True
        result.note = ("RetroRules suggestion (rule-based prediction — verify each step): "
                       + result.note)
        # For rule-based reactions the suggested id defaults to the reaction NAME (the
        # user can still edit it in the Apply dialog) — the auto id from the SMILES-slug
        # id would be meaningless.
        if not result.reactions.empty and "name" in result.reactions.columns:
            for i in result.reactions.index:
                nm = str(result.reactions.at[i, "name"] or "").strip()
                if nm:
                    result.reactions.at[i, "suggested_id"] = nm
        self.pathway_panel.append_result(result)
        self._pathway_results = getattr(self, "_pathway_results", [])
        self._pathway_results.append(result)
        QMessageBox.information(
            self, "Added as a suggested pathway",
            "The rule-based steps are now a suggested pathway in Pathway Design. Review "
            "them, then click “Apply pathway” to add them to your model — you can rename "
            "the reactions and metabolites in the dialog that opens.\n\n"
            "Remember these are predictions to verify, not database reactions.")

    def _predict_pathway(self, target: str, min_flux: float,
                         forbidden: set | None = None) -> None:
        """Search for routes to `target`. `forbidden` excludes reactions from the
        search — used by "Find a route that runs" to steer around a blocking step."""
        if self.project is None:
            return
        if not self._has_loaded_db():
            QMessageBox.information(
                self, "No reaction database loaded",
                "Pathway prediction needs a reaction database. In the “Reaction databases” "
                "panel, tick one or more databases and click “Load selected databases”, then "
                "run Predict pathway again.")
            return
        from ..core.network_graph import short_metabolite_name
        from .dialogs.pathway_settings_dialog import PathwaySettingsDialog

        # Offer the start-metabolite picker from the host model's metabolites.
        met_items = [(short_metabolite_name(m.id, m.name or "") + f"  ({m.id})", m.id)
                     for m in self.project.model.metabolites]
        met_items.sort(key=lambda t: t[0].lower())
        target_display = self.pathway_panel.target_combo.currentText() or target
        dlg = PathwaySettingsDialog(self, target_display, met_items, min_flux=min_flux)
        if dlg.exec() != PathwaySettingsDialog.Accepted:
            return
        cfg = dlg.values()

        host = self.project.model.copy()
        start = None if cfg["use_all"] else (cfg["start_metabolites"] or None)
        n_alt = cfg["n_alternatives"]
        pref_ec = cfg["preferred_ec"] or None
        max_steps = cfg["max_steps"]
        algorithm = cfg.get("algorithm", "retro")
        priority = cfg.get("priority", "yield")
        include_boundary = cfg.get("include_boundary", False)
        flux_only = cfg.get("flux_carrying_starts_only", False)

        # Restricting the starting pool to compounds the host actually makes is all it
        # takes to stop a route being anchored on an idle metabolite — the search then has
        # to find the heterologous steps that supply it. Done here rather than inside the
        # search so the default behaviour is untouched for anyone who leaves it off.
        flux_note = ""
        if flux_only:
            from ..core import flux_context as fx
            starts, ctx = fx.flux_carrying_starts(host, requested=start)
            if not ctx.usable:
                QMessageBox.warning(
                    self, "Could not identify flux-carrying precursors",
                    ctx.summary() + "<br><br>Searching without the restriction instead.")
            else:
                start = starts
                flux_note = ctx.summary()

        def work(report_progress=None):
            # Build the combined universal now (merge is cached in self._pathway_db and
            # only rebuilt when the selection changes), resolve the target, then search.
            combined = self._combined_included_model()
            resolved = self._resolve_target_in(combined, target)
            def _progress(done, total, message):
                # Report which alternative is being searched, so a long multi-route
                # search is not a silent wait (VI.14b).
                if report_progress is not None:
                    report_progress(f"{message}  ({done}/{total})")

            results = pathway_design.find_pathways(
                host, resolved, combined, start_metabolites=start, max_steps=max_steps,
                n_alternatives=n_alt, preferred_ec=pref_ec, algorithm=algorithm,
                priority=priority, include_boundary=include_boundary,
                forbidden_reactions=forbidden or None, progress=_progress)
            # Issue 14: namespace-match diagnostics restricted to the route metabolites.
            report = None
            try:
                from ..core import namespace
                relevant = set()
                for r in results:
                    for rid in r.reaction_ids:
                        if combined.reactions.has_id(rid):
                            relevant |= {m.id for m in combined.reactions.get_by_id(rid).metabolites}
                if relevant:
                    report = namespace.translation_report(host, combined, relevant_ids=relevant)
            except Exception:  # noqa: BLE001
                report = None
            return results, report

        def finish(payload):
            results, report = payload
            if flux_note:
                # Say the search was restricted, and to what — otherwise a route with an
                # unexpected extra step (the mgsA-type bridge) looks arbitrary.
                for r in results:
                    r.note = (r.note + "  Starting points were restricted to "
                              "flux-carrying metabolites. " + flux_note).strip()
            self._pathway_results = results
            self.pathway_panel.show_results(results)
            found = sum(1 for r in results if not r.reactions.empty)
            self.status_label.setText(
                f"Pathway prediction complete for {target}: {found} pathway(s) found.")
            self._show_namespace_report(report)
            # "Not found" is the least useful answer we can give: very often a close
            # relative IS reachable and the real target is one step beyond it (L6).
            if not found:
                self._offer_reachable_analogues(target, combined)

        ok, res = run_busy(
            self, f"Predicting pathway to {target}…", work,
            title="Pathway design", after=finish,
            after_message="Preparing results…", cancelable=True, progress=True)
        if not ok and not was_cancelled(res):
            QMessageBox.critical(self, "Pathway prediction failed", str(res))

    def _explore_upstream(self) -> None:
        """What feeds this route's entry into native metabolism, and what else could.

        The hand-over point decides whether a design runs at all — an entry compound that
        carries no flux (lactaldehyde in Synechocystis) makes the route unbuildable until
        something supplies it. Purely additive: it reads the finished route.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        ids = list(result.reaction_ids)
        host = self.project.model

        def work():
            from ..core import upstream as up
            return up.explore_upstream(host, db, ids)

        ok, report = run_busy(
            self, "Looking at what feeds this route…", work,
            title="Upstream exploration", cancelable=True)
        if not ok:
            if not was_cancelled(report):
                QMessageBox.warning(self, "Could not explore upstream", str(report))
            return
        if not report.entry_metabolites:
            QMessageBox.information(
                self, "No entry point",
                "This route does not draw any non-cofactor compound from your host, so "
                "there is no hand-over point to explore upstream of.")
            return
        if not report.total_candidates:
            QMessageBox.information(
                self, "Nothing upstream found",
                report.headline().replace("<b>", "").replace("</b>", "")
                + "\n\nNo database reaction produces those compounds. Load more "
                  "databases, or fetch chemistry around them.")
            return

        from .dialogs.upstream_dialog import UpstreamDialog
        dlg = UpstreamDialog(self, report)
        dlg.reactions_chosen.connect(
            lambda picked, r=result: self._extend_route_upstream(r, picked))
        dlg.exec()

    def _extend_route_upstream(self, result, reaction_ids) -> None:
        """Add chosen upstream reactions to a route, as a new tab.

        A copy rather than an edit in place: the original design stays intact for
        comparison, which is the whole point of proposing an alternative supply.
        """
        db = self._combined_included_model()
        if db is None or not reaction_ids:
            return
        extended = result.duplicate()
        added = [rid for rid in reaction_ids if rid not in extended.reaction_ids]
        if not added:
            return
        extended.reaction_ids = list(extended.reaction_ids) + added
        try:
            rebuilt = pathway_design._result_from_ids(
                self.project.model, db, extended.target, extended.reaction_ids)
            rebuilt.note = (f"Extended upstream with {len(added)} reaction(s): "
                            + ", ".join(added) + ". " + rebuilt.note)
            self.pathway_panel.append_result(rebuilt)
            self.status_label.setText(
                f"Added {len(added)} upstream reaction(s) to the route.")
        except Exception as exc:  # noqa: BLE001 - report rather than lose the selection
            QMessageBox.warning(self, "Could not extend the route", str(exc))

    def _structural_search(self) -> None:
        """Find a target by structure when its exact name is not in the database (VI.2)."""
        db = self._combined_included_model()
        if db is None:
            QMessageBox.information(
                self, "No database loaded",
                "Load a reaction database first — structural search looks for compounds "
                "inside the databases you have loaded.")
            return
        from .dialogs.structural_search_dialog import StructuralSearchDialog
        dlg = StructuralSearchDialog(self, db)
        dlg.target_chosen.connect(self._select_target_id)
        dlg.exec()

    def _select_target_id(self, met_id: str) -> None:
        """Point the target selector at a metabolite id chosen elsewhere."""
        combo = getattr(self.pathway_panel, "target_combo", None)
        if combo is None or not met_id:
            return
        for i in range(combo.count()):
            text = combo.itemText(i)
            if met_id == combo.itemData(i) or f"({met_id})" in text:
                combo.setCurrentIndex(i)
                self.status_label.setText(f"Target set to {met_id}.")
                return
        # Not in the list (it may live only in the database, not the host): add it.
        combo.addItem(met_id, met_id)
        combo.setCurrentIndex(combo.count() - 1)
        self.status_label.setText(f"Target set to {met_id}.")

    def _feasibility_report(self) -> None:
        """Everything known about whether the selected route is buildable (VI.10).

        Runs the checks that have not been run yet — isomer check with structures fetched
        on demand, blockage analysis — and offers one-click gap filling when a blocking
        metabolite is identified.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        ids = result.reaction_ids

        def build():
            from ..core import chemistry, feasibility as fz, gapfill
            from ..core import pathway_diagnostics as pdg
            eng = pathway_design.apply_pathway(self.project.model, db, ids,
                                               rename=False, target=result.target)
            native = {m.id for m in self.project.model.metabolites}
            # Fetch whatever chemistry is missing FIRST, rather than reporting "cannot be
            # checked" and leaving the user to go and find it: structures come from each
            # metabolite's cross-references, and a formula is derived from the structure
            # where the database has none. Both checks below then have something to work
            # with, and the result is cached so this costs nothing next time.
            backfill = chemistry.resolve_missing_chemistry(db, ids, online=True)
            if backfill.formulas_added:
                # Balance could not be verified before; with formulas it can be.
                self._recheck_balance(result, db, ids)
            # Re-run the isomer check WITH structure fetching, which the search itself
            # deliberately skips to stay fast.
            iso = chemistry.check_pathway_isomers(db, ids, online=True)
            result.isomer_warnings = [f.as_sentence() for f in iso.findings]
            result.isomer_checked = bool(iso.checked)
            result.isomer_coverage = iso.coverage
            try:
                diag = pdg.diagnose(eng, ids, target_id=result.target, native_ids=native)
            except Exception:  # noqa: BLE001
                diag = None
            blockage = gapfill.find_blockers(eng, ids, native_ids=native)
            if blockage.blocked and diag is not None and blockage.blockers:
                # Prefer the operational blocker list over the topological one (VI.14c).
                diag.blocking_metabolites = blockage.blocker_labels
                diag.recommendation = blockage.recommendation
            rep = fz.assess(result, diagnosis=diag,
                            branching=getattr(result, "_branching", None),
                            mdf=getattr(result, "_mdf", None))
            rep.backfill_note = backfill.sentence()
            return rep, blockage

        ok, built = run_busy(self, "Assessing this route…", build,
                             title="Feasibility", cancelable=True)
        if not ok:
            if not was_cancelled(built):
                QMessageBox.warning(self, "Could not assess", str(built))
            return
        rep, blockage = built
        self.pathway_panel.refresh_action_visibility()

        from .dialogs.feasibility_dialog import FeasibilityDialog
        FeasibilityDialog(
            self, rep,
            blockage=blockage,
            on_fetch_missing=(lambda b=blockage: self._autofill_gap(b))
            if blockage.blockers else None).exec()

    @staticmethod
    def _route_metabolites(db, reaction_ids):
        """``(id, name)`` for the compounds whose concentration the MDF can set.

        Water and protons are excluded: their activity is folded into ΔrG′° rather than
        being a free variable, so offering to set them would be misleading.
        """
        out, seen = [], set()
        for rid in reaction_ids:
            if not db.reactions.has_id(rid):
                continue
            for met in db.reactions.get_by_id(rid).metabolites:
                if met.id in seen:
                    continue
                seen.add(met.id)
                base = (met.name or met.id).lower()
                if base.startswith(("h2o", "water")) or base in ("h+", "proton") \
                        or met.id.lower().startswith(("h2o", "h_")):
                    continue
                out.append((met.id, met.name or met.id))
        return sorted(out, key=lambda p: p[1].lower())

    @staticmethod
    def _recheck_balance(result, db, ids) -> None:
        """Re-run the balance check after formulas were back-filled.

        Steps recorded as "could not be verified" may now be verifiable — and some of them
        may turn out to be genuinely unbalanced. Both outcomes are more useful than the
        original "unknown", so the result is updated in place.
        """
        unverified, reasons, unbalanced = [], {}, False
        for rid in ids:
            if not db.reactions.has_id(rid):
                continue
            verdict = pathway_design.reaction_balance_verdict(db.reactions.get_by_id(rid))
            if not verdict.checkable:
                unverified.append(rid)
                reasons[rid] = list(verdict.missing_formula_labels)
            elif verdict.verifiably_imbalanced:
                unbalanced = True
        result.unverified_steps = unverified
        result.unverified_reasons = reasons
        result.balanced = not unbalanced

    def _autofill_gap(self, blockage) -> None:
        """Fetch the chemistry a blocked route is missing and complete it in one step.

        Fetched reactions are used for this session only; the user is offered the chance
        to save them deliberately rather than having their merged database mutated.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None:
            return
        names = ", ".join(b.label for b in blockage.blockers[:3])
        if QMessageBox.question(
                self, "Fetch missing reactions?",
                f"This route is blocked at <b>{blockage.bottleneck}</b> because the host "
                f"cannot supply <b>{names}</b>.<br><br>"
                "Fetch reactions around it from KEGG and Rhea, then search again for a "
                "complete route?<br><br>"
                "<span style='color:#5f6368'>This can add several hundred reactions and "
                "may take a minute. They are used for this session only.</span>",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes) != QMessageBox.Yes:
            return

        from ..core import gapfill

        def work(report=None):
            return gapfill.autofill(self.project.model, db, result.target,
                                    blockage.blockers)

        ok, res = run_busy(self, "Fetching the missing chemistry…", work,
                           title="Fetch missing reactions", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                QMessageBox.warning(self, "Could not fetch", str(res))
            return
        if not res.success:
            QMessageBox.information(self, "No complete route found", res.message)
            return

        if QMessageBox.question(
                self, "Complete route found",
                f"{res.message}<br><br>Show this completed route as the current result?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes) != QMessageBox.Yes:
            return
        try:
            new_res = pathway_design._result_from_ids(
                self.project.model, res.merged_model, result.target,
                res.reaction_ids, rank=0)
            new_res.note = ("Completed by fetching missing chemistry: " + new_res.note)
            self._pathway_db = res.merged_model
            self.pathway_panel.append_result(new_res)
            self._pathway_results = getattr(self, "_pathway_results", [])
            self._pathway_results.append(new_res)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Could not show the route", str(exc))
            return
        QMessageBox.information(
            self, "Added",
            "The completed route is now the current result.<br><br>"
            "The fetched reactions are loaded for this session. Use <b>Save loaded "
            "database…</b> in the Reaction databases panel if you want to keep them.")

    def _declare_final_step(self) -> None:
        """Let the user add ONE declared enzyme step from this route's product to their
        real target, then re-evaluate the completed route (L6).

        The database frequently stops one well-known reaction short of the goal —
        N-methyltryptamine is present but DMT is not, though the same methyltransferase
        performs both methylations. Rather than reporting failure, the user declares that
        step (naming the enzyme and any cofactors) and the whole pathway becomes
        analysable end-to-end.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        src = result.target
        name, ok = QInputDialog.getText(
            self, "Add a final step",
            f"This route makes <b>{src}</b>.\n\n"
            "Name the compound your final enzyme converts it into\n"
            "(e.g. “N,N-dimethyltryptamine”):")
        if not ok or not name.strip():
            return
        product = name.strip()
        enzyme, ok = QInputDialog.getText(
            self, "Add a final step",
            f"Which enzyme converts {src} → {product}?\n"
            "(name or EC number; used to label the reaction)")
        if not ok:
            return
        enzyme = enzyme.strip() or "declared final step"
        cof, ok = QInputDialog.getText(
            self, "Add a final step",
            "Optional co-substrate/co-product pairs consumed by this step,\n"
            "as model metabolite ids separated by commas\n"
            "(e.g. “amet_c>ahcys_c” for a SAM-dependent methylation). Leave blank if none:")
        if not ok:
            cof = ""

        def build():
            import cobra
            from ..core import pathway_design as pdz
            eng = pdz.apply_pathway(self.project.model, db, result.reaction_ids,
                                    rename=False, target=src)
            if not eng.metabolites.has_id(src):
                raise RuntimeError(f"'{src}' is not present after applying the route.")
            src_met = eng.metabolites.get_by_id(src)
            pid = "".join(c if c.isalnum() else "_" for c in product).strip("_")[:40]
            pid = f"{pid}_c"
            if eng.metabolites.has_id(pid):
                prod_met = eng.metabolites.get_by_id(pid)
            else:
                prod_met = cobra.Metabolite(pid, name=product, compartment="c")
                eng.add_metabolites([prod_met])
            rid = "DECLARED_" + pid[:-2].lower()
            rxn = cobra.Reaction(rid, name=f"{product} synthesis ({enzyme})")
            rxn.lower_bound, rxn.upper_bound = 0.0, 1000.0
            coeffs = {src_met: -1.0, prod_met: 1.0}
            for pair in (cof or "").split(","):
                pair = pair.strip()
                if not pair:
                    continue
                a, _, b = pair.partition(">")
                for mid, sign in ((a.strip(), -1.0), (b.strip(), 1.0)):
                    if mid and eng.metabolites.has_id(mid):
                        m = eng.metabolites.get_by_id(mid)
                        coeffs[m] = coeffs.get(m, 0.0) + sign
            rxn.add_metabolites(coeffs)
            rxn.annotation["gsm.declared"] = "user-declared final step"
            if enzyme:
                rxn.annotation["ec-code"] = enzyme
            eng.add_reactions([rxn])
            sink = self._resolve_product_sink_impl(eng, pid)
            flux = float("nan")
            if sink:
                eng.objective = sink
                v = eng.slim_optimize()
                flux = float(v) if v is not None and v == v else float("nan")
            return eng, pid, rid, flux

        ok2, built = run_busy(self, f"Completing the route to {product}…", build,
                              title="Add a final step", cancelable=True)
        if not ok2:
            if not was_cancelled(built):
                QMessageBox.warning(self, "Could not add the step", str(built))
            return
        _eng, pid, rid, flux = built
        carries = flux == flux and flux > 1e-9
        QMessageBox.information(
            self, "Final step added",
            f"Declared <b>{rid}</b>: {src} → {product} ({enzyme}).<br><br>"
            + (f"The completed route <b>carries {flux:.4g}</b> mmol gDW⁻¹ h⁻¹ to "
               f"{product}." if carries else
               "The completed route <b>carries no flux</b> — the limitation is upstream "
               "of this final step, not in it.")
            + "<br><br><span style='color:#5f6368'>This step is a declaration, not a "
              "database reaction: you are asserting the enzyme exists and accepts this "
              "substrate. Verify that before relying on the result.</span>")

    def _offer_reachable_analogues(self, target: str, db) -> None:
        """After a failed search, suggest structurally/chemically related compounds that
        the host CAN reach, and offer to search for one of them instead (L6).

        This is how "violacein is not in the database" becomes "violaceinate is, with 9
        producing reactions", and how a dead end turns into a design that is one known
        enzyme short of the goal.
        """
        from ..core import pathway_search
        name = self._plain_compound_name(target)
        smiles = ""
        try:
            met = None
            for mdl in (self.project.model if self.project else None, db):
                if mdl is not None and mdl.metabolites.has_id(target):
                    met = mdl.metabolites.get_by_id(target)
                    break
            if met is not None:
                from ..core import retrorules
                smiles = retrorules.metabolite_smiles(met) or ""
        except Exception:  # noqa: BLE001
            smiles = ""

        ok, hits = run_busy(
            self, "Looking for related compounds the host can reach…",
            lambda: pathway_search.nearest_reachable_analogues(
                db, name, target_smiles=smiles, limit=6),
            title="No route found", cancelable=True)
        if not ok or not hits:
            return
        lines = "".join(
            f"<li><b>{h['name']}</b> <span style='color:#5f6368'>({h['id']}, "
            f"{h['n_producers']} producing reaction(s)) — {h['reason']}</span></li>"
            for h in hits)
        box = QMessageBox(self)
        box.setWindowTitle("No route found — but these are reachable")
        box.setIcon(QMessageBox.Information)
        box.setTextFormat(Qt.RichText)
        box.setText(
            f"No route to <b>{name}</b> was found, but the database contains related "
            f"compounds that <i>can</i> be produced:<ul>{lines}</ul>"
            "Searching for one of these often gets you within a single known enzyme of "
            "your target.<br><br><span style='color:#5f6368'>These are suggested by "
            "structure and chemical name similarity — check that the chemistry is "
            "genuinely related before relying on one.</span>")
        search_btn = box.addButton("Search for the first suggestion",
                                   QMessageBox.AcceptRole)
        box.addButton("Close", QMessageBox.RejectRole)
        box.exec()
        if box.clickedButton() is search_btn:
            self._predict_pathway(hits[0]["id"], self.pathway_panel.min_flux.value())

    def _apply_pathway(self) -> None:
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None:
            return
        ids = result.reaction_ids
        if not ids:
            return
        # User-edited "suggested id"s become the reaction ids on apply (#B5).
        rename_map = {}
        if not result.reactions.empty and {"reaction", "suggested_id"} <= set(result.reactions.columns):
            rename_map = {str(r): str(s) for r, s in
                          zip(result.reactions["reaction"], result.reactions["suggested_id"])
                          if str(s).strip()}

        # Dry-run on a copy to summarise exactly what will change (#B9).
        model = self.project.model
        before_mets = {m.id for m in model.metabolites}
        rep: dict = {}
        try:
            engineered = pathway_design.apply_pathway(
                model, db, ids, rename=True, rename_map=rename_map, report=rep,
                target=result.target)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Could not apply pathway", str(exc))
            return
        added_ids = rep.get("added_reactions", [])
        rxn_rows = [{"id": rid, "name": engineered.reactions.get_by_id(rid).name or "",
                     "equation": engineered.reactions.get_by_id(rid).build_reaction_string(
                        use_metabolite_names=True)}
                    for rid in added_ids if engineered.reactions.has_id(rid)]
        new_mets = [{"id": m.id, "name": m.name}
                    for m in engineered.metabolites if m.id not in before_mets]
        target_met = pathway_search._resolve_target(engineered, engineered, result.target)
        can_secrete = any(c.startswith("e") for c in (getattr(engineered, "compartments", {}) or {}))

        from .dialogs.add_pathway_dialog import AddPathwayDialog
        dlg = AddPathwayDialog(
            self, target_name=self._pathway_target_label(result.target),
            reactions=rxn_rows, new_metabolites=new_mets,
            categories=self.project.categories.names(),
            collapsed=rep.get("collapsed_compartments", []),
            dropped=rep.get("dropped_reactions", []), unbalanced=rep.get("unbalanced", []),
            unverified=rep.get("unverified", []),
            can_secrete=can_secrete and target_met is not None,
            product_exchange=rep.get("product_exchange"))
        if dlg.exec() != AddPathwayDialog.Accepted:
            return
        opts = dlg.values()

        before = {r.id for r in model.reactions}
        rep2: dict = {}
        route_id_holder = {}

        def _do(m):
            # add_pathway_to_model guarantees a single product route (#F3); the route
            # kind (exported / gas / intracellular / none) is the user's choice (#B).
            pathway_design.add_pathway_to_model(
                m, db, ids, rename=True, rename_map=rename_map, report=rep2,
                target=result.target, product_route=opts.get("route", "secrete"))
            px = rep2.get("product_exchange")
            if px and px.get("id"):
                route_id_holder["id"] = px["id"]
            # Apply the user's manual id/name edits to the freshly-added objects (#B).
            self._apply_pathway_edits(m, opts)

        ok, err = run_busy(
            self, f"Adding {len(ids)} heterologous reaction(s) to the model…",
            lambda: self.project.apply_edit(_do), title="Applying pathway", cancelable=True)
        if not ok:
            if not was_cancelled(err):
                QMessageBox.warning(self, "Could not apply pathway", str(err))
            return
        added = sorted({r.id for r in model.reactions} - before)
        try:
            result.applied_ids = added
        except Exception:  # noqa: BLE001
            self._applied_by_target[result.target] = added
        # Assign to the chosen category.
        cat = opts.get("category", "").strip()
        if cat and added:
            if not self.project.categories.has(cat):
                self.project.categories.create(cat)
            self._add_ids_to_category(cat, added)
        self._refresh_after_edit(
            f"Applied pathway: added {len(added)} reaction(s)"
            + (f" + export route {route_id_holder['id']}" if route_id_holder.get("id") else "")
            + (f" · category '{cat}'" if cat else "") + ".")

    @staticmethod
    def _apply_pathway_edits(m, opts: dict) -> None:
        """Apply the user's manual id/name edits (from the add-pathway dialog) to the
        freshly-added reactions and metabolites. Ids are set via the private ``_id``
        with a single ``model.repair()`` at the end (avoids O(n^2) reindexing), and
        never overwrite an existing id."""
        changed = False
        for kind, coll in (("reactions", m.reactions), ("metabolites", m.metabolites)):
            renames = opts.get(f"{'rxn' if kind == 'reactions' else 'met'}_renames", {})
            names = opts.get(f"{'rxn' if kind == 'reactions' else 'met'}_names", {})
            for orig_id in set(renames) | set(names):
                if not coll.has_id(orig_id):
                    continue
                obj = coll.get_by_id(orig_id)
                new_name = names.get(orig_id)
                if new_name and new_name != (obj.name or ""):
                    obj.name = new_name
                new_id = renames.get(orig_id)
                if new_id and not coll.has_id(new_id):
                    obj._id = new_id
                    changed = True
        if changed:
            m.repair()

    # ----- growth settings (Issue 5) -----------------------------------
    def _apply_growth_settings(self, medium: dict, biomass_id) -> None:
        if self.project is None:
            return

        def _edit(m):
            if biomass_id and m.reactions.has_id(biomass_id):
                m.objective = biomass_id
            if medium:
                # Keep only exchange ids that exist; set as the medium.
                valid = {k: v for k, v in medium.items() if m.reactions.has_id(k)}
                m.medium = valid
        self.project.apply_edit(_edit)
        self._refresh_after_edit("Applied growth settings.")

    def _apply_growth_mode(self, mode: str) -> None:
        if self.project is None:
            return
        from ..core import physiology
        presets = {
            "autotrophic": physiology.apply_photoautotrophic_preset,
            "mixotrophic": physiology.apply_mixotrophic_preset,
            "heterotrophic": physiology.apply_heterotrophic_mode,
        }
        fn = presets.get(mode)
        if fn is None:
            return
        notes = self.project.apply_edit(lambda m: fn(m))
        self._refresh_after_edit(f"Applied {mode} growth mode.")

        if mode == "autotrophic":
            QMessageBox.information(self, "Autotrophic growth mode",
                                    "\n".join(f"• {n}" for n in notes))
            return

        # Mixotrophy and heterotrophy are *defined* by having organic carbon, but which
        # carbon is an experimental choice the preset cannot make — a model may offer
        # thirty candidates. Ask now, with the options listed, instead of reporting what
        # was done and leaving the user to find the right exchange in the medium editor.
        from .dialogs.carbon_source_dialog import CarbonSourceDialog

        chosen = CarbonSourceDialog.choose(self, self.project.model, mode, notes)
        if not chosen:
            return

        def _feed(model):
            opened = []
            for exchange, rate in chosen.items():
                if model.reactions.has_id(exchange):
                    model.reactions.get_by_id(exchange).lower_bound = -abs(rate)
                    opened.append(f"{exchange} at {rate:g} mmol gDW⁻¹ h⁻¹")
            return opened

        opened = self.project.apply_edit(_feed)
        self._refresh_after_edit(
            f"Applied {mode} growth mode with {len(opened)} carbon source(s).")
        QMessageBox.information(
            self, f"{mode.capitalize()} growth mode",
            "\n".join(f"• {n}" for n in list(notes) + opened))

    # ----- Strategy Explorer (graphical engine) ------------------------
    def _save_strategy(self) -> None:
        if self.project is None:
            return
        result = getattr(self, "_last_result", None)
        if result is None or not getattr(result, "is_optimal", False) or result.fluxes.empty:
            QMessageBox.information(
                self, "Save strategy",
                "Run FBA or pFBA first — a strategy captures the current model together "
                "with its solved fluxes, so a flux state is needed before saving.")
            return
        n = len(self.project.strategies) + 1
        name, ok = QInputDialog.getText(self, "Save flux state as strategy",
                                        "Name this strategy (round of engineering):",
                                        text=f"Round {n}")
        if not ok or not name.strip():
            return
        # Optional product reaction for the titre waterfall. Offer EVERY reaction, not
        # just exchanges: a heterologous product added by Pathway Design often has no
        # exchange of its own, so an exchange-only list silently omitted the very
        # compound the user was engineering. Exchanges stay first (the common case),
        # and the field is editable so an id can simply be typed.
        # `.exchanges` raises on models whose external compartment cannot be guessed;
        # never let that stop a strategy being saved.
        try:
            ex_ids = sorted(r.id for r in self.project.model.exchanges)
        except Exception:  # noqa: BLE001
            ex_ids = sorted(r.id for r in self.project.model.reactions if r.boundary)
        other_ids = sorted(r.id for r in self.project.model.reactions
                           if r.id not in set(ex_ids))
        choices = ["(none)"] + ex_ids + other_ids
        target, ok = QInputDialog.getItem(
            self, "Reaction of interest (optional)",
            "Reaction to follow in the titre waterfall (optional) — exchanges first, "
            "then all other reactions; you can also type an id:",
            choices, 0, True)
        target = "" if (not ok or target == "(none)") else target.strip()
        from ..core.flux_state import FluxState
        state = FluxState(
            name=name.strip(), fluxes={k: float(v) for k, v in result.fluxes.items()},
            objective_value=float(result.objective_value), method=result.method,
            target=target)
        self.project.strategies.add(state)
        self.project._dirty = True
        self.strategy_explorer.set_strategies(self.project.strategies)
        self.escher_explorer.set_strategies(self.project.strategies)
        self.status_label.setText(f"Saved strategy '{state.name}'.")

    def _remove_strategy(self, name: str) -> None:
        if self.project is None:
            return
        self.project.strategies.remove(name)
        self.project._dirty = True
        self.strategy_explorer.set_strategies(self.project.strategies)
        self.escher_explorer.set_strategies(self.project.strategies)

    def _open_escher(self) -> None:
        """Open the interactive Escher map tab and build the map for the current model.

        The explorer is added to the tab strip on demand — it is reached from the Network
        Visualization menu rather than living there permanently. Without the insert,
        `setCurrentWidget` is a silent no-op on a widget the tab bar does not own, so the
        map was built (a few seconds of frozen UI) into a widget nobody could see.
        """
        if self.tabs.indexOf(self.escher_explorer) < 0:
            self.tabs.addTab(self.escher_explorer, "Escher Visualizer")
        self.tabs.setCurrentWidget(self.escher_explorer)
        if self.project is None:
            return
        # Building a map is seconds of synchronous work; say so rather than looking hung.
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self.escher_explorer.rebuild()
        finally:
            QApplication.restoreOverrideCursor()

    def _open_plot_gallery(self) -> None:
        """Assemble the publication-grade figures applicable to the current state
        and open the Plot Gallery (SVG/PDF/CSV export)."""
        if self.project is None:
            return
        from .dialogs.plot_gallery_dialog import PlotGalleryDialog
        from .viz import plots
        df = getattr(self, "_last_table", None)
        cols = set(df.columns) if df is not None and not getattr(df, "empty", True) else set()
        entries = []

        if {"minimum", "maximum"} & cols and {"minimum", "maximum"}.issubset(cols):
            entries.append(("FVA flux ranges (tornado)",
                            lambda c: c.render(plots.fva_tornado, df, title="fva_ranges")))
        if {"target_type", "flux_start", "flux_end"}.issubset(cols):
            entries.append(("FSEOF scan",
                            lambda c: c.render(plots.fseof_scan, df, title="fseof_scan")))
        if {"predicted_growth"} & cols or {"guaranteed_product"} & cols:
            entries.append(("Strain-design comparison",
                            lambda c: c.render(plots.strain_design_comparison, df,
                                               title="design_comparison")))
        # Exchange-flux bars from the last solved FBA/pFBA state.
        result = getattr(self, "_last_result", None)
        if result is not None and getattr(result, "is_optimal", False) and not result.fluxes.empty:
            fluxes = {k: float(v) for k, v in result.fluxes.items()}
            ex_ids = [r.id for r in self.project.model.exchanges]
            labels = {r.id: (r.name or r.id) for r in self.project.model.exchanges}
            entries.append(("Exchange fluxes (uptake/secretion)",
                            lambda c: c.render(plots.exchange_flux_bars, fluxes, ex_ids, labels,
                                               title="exchange_fluxes")))
        # Strategy figures, if any strategies are saved.
        strat = self.project.strategies
        if len(strat) >= 1:
            names = strat.names()
            ids = sorted({k for s in strat for k, v in s.fluxes.items() if abs(v) > 1e-6})[:25]
            matrix = strat.flux_matrix(ids)
            entries.append(("Multi-strategy heatmap",
                            lambda c: c.render(plots.multi_strategy_heatmap, matrix, names,
                                               title="strategy_heatmap")))
            wf = [(s.target_flux() if s.target_flux() == s.target_flux() else 0.0) for s in strat]
            entries.append(("Titre waterfall",
                            lambda c: c.render(plots.titre_waterfall, names, wf,
                                               title="titre_waterfall")))

        if not entries:
            QMessageBox.information(
                self, "Plot Gallery",
                "No figures are available for the current result yet. Run an analysis "
                "(FBA/pFBA, FVA, FSEOF, strain design…) or save some strategies, then "
                "open the Plot Gallery.")
            return
        PlotGalleryDialog(self, entries).exec()

    def _show_namespace_report(self, report) -> None:
        """Issue 14: warn when route metabolites matched the host ambiguously or
        several database ids collapsed onto one host metabolite."""
        if not report:
            return
        ambiguous = report.get("ambiguous") or []
        collisions = report.get("collisions") or []
        if not (ambiguous or collisions):
            return
        lines = []
        if ambiguous:
            detail = "; ".join(f"{db} → {', '.join(hosts)}" for db, hosts in ambiguous[:6])
            lines.append("Ambiguous matches (a database metabolite whose cross-references "
                         f"point at more than one host metabolite): {detail}"
                         + (" …" if len(ambiguous) > 6 else ""))
        if collisions:
            detail = "; ".join(f"{host} ← {', '.join(dbs)}" for host, dbs in collisions[:6])
            lines.append("Collisions (several database metabolites mapped onto the same "
                         f"host metabolite): {detail}"
                         + (" …" if len(collisions) > 6 else ""))
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("Namespace match — please review")
        box.setText("The pathway was attached to your model's metabolites, but some "
                    "matches were not unambiguous:")
        box.setInformativeText("\n\n".join(lines)
                               + "\n\nConfirm the pathway connects to the intended "
                               "metabolites before relying on it.")
        box.exec()

    def _show_pathway_apply_report(self, rep: dict) -> None:
        """Issue 7: after applying a pathway, warn about compartment collapses,
        dropped (pure-transport) reactions and any mass/charge-unbalanced steps."""
        if not rep:
            return
        collapsed = rep.get("collapsed_compartments") or []
        dropped = rep.get("dropped_reactions") or []
        unbalanced = rep.get("unbalanced") or []
        unverified = rep.get("unverified") or []
        if not (collapsed or dropped or unbalanced or unverified):
            return
        lines = []
        if collapsed:
            lines.append("Compartments collapsed into the host's main compartment: "
                         + ", ".join(collapsed)
                         + "\n(metabolites from compartments your host lacks were merged — "
                         "review that no chemically distinct species were unified).")
        if dropped:
            lines.append(f"{len(dropped)} reaction(s) became pure transport after the "
                         f"collapse and were dropped: {', '.join(dropped[:8])}"
                         + (" …" if len(dropped) > 8 else ""))
        if unbalanced:
            detail = "; ".join(f"{rid} ({res})" for rid, res in unbalanced[:8])
            lines.append(f"⚠ {len(unbalanced)} added reaction(s) are mass/charge-"
                         f"unbalanced: {detail}"
                         + (" …" if len(unbalanced) > 8 else "")
                         + "\nUse Model QC / the reaction editor to balance them.")
        if unverified:
            # Deliberately NOT counted with the unbalanced ones: nothing is wrong with
            # these reactions — a participant simply has no formula, so their atoms were
            # never counted. Reporting them together taught users to ignore both.
            detail = "; ".join(f"{rid} ({why})" for rid, why in unverified[:6])
            lines.append(f"{len(unverified)} added reaction(s) could not be balance-"
                         f"checked: {detail}"
                         + (" …" if len(unverified) > 6 else "")
                         + "\nThis is a gap in the database, not an imbalance. The "
                         "reaction editor can fetch the missing formula.")
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning if unbalanced else QMessageBox.Information)
        box.setWindowTitle("Pathway applied — balance report")
        box.setText("The pathway was applied, with the following notes:")
        box.setInformativeText("\n\n".join(lines))
        box.exec()

    def _pathway_target_label(self, target: str) -> str:
        """A readable target name for a category default (e.g. 'Astaxanthin')."""
        m = (self.project.model.metabolites.get_by_id(target)
             if self.project and self.project.model.metabolites.has_id(target) else None)
        if m is not None and m.name:
            return m.name
        from ..core.network_graph import short_metabolite_name
        label = short_metabolite_name(target, "")
        return label.rsplit("_", 1)[0] if "_" in label else label

    def _prompt_pathway_category(self, target: str, ids: list) -> None:
        if not ids or self.project is None:
            return
        default = f"{self._pathway_target_label(target)} production"
        existing = [n for n in self.project.categories.names() if n != default]
        name, ok = QInputDialog.getItem(
            self, "Group the added reactions",
            "Add the added reaction(s) to a category (pick an existing one or type a new name):",
            [default] + existing, 0, True)
        if not ok or not name.strip():
            return
        name = name.strip()
        if not self.project.categories.has(name):
            self.project.categories.create(name)
        self._add_ids_to_category(name, ids)

    def _applied_ids_for(self, result) -> list:
        return list(getattr(result, "applied_ids", None)
                    or self._applied_by_target.get(getattr(result, "target", ""), []))

    def _remove_pathway(self) -> None:
        result = self.pathway_panel.current_result()
        if self.project is None or result is None:
            return
        ids = [rid for rid in self._applied_ids_for(result)
               if self.project.model.reactions.has_id(rid)]
        if not ids:
            QMessageBox.information(
                self, "Remove pathway",
                "This pathway hasn't been added to the model (nothing to remove).")
            return
        if QMessageBox.question(
                self, "Remove pathway from model",
                f"Remove the {len(ids)} reaction(s) added for this pathway?") != QMessageBox.Yes:
            return

        def _remove(m):
            m.remove_reactions([m.reactions.get_by_id(r) for r in ids
                                if m.reactions.has_id(r)], remove_orphans=True)
        self.project.apply_edit(_remove)
        if hasattr(result, "applied_ids"):
            result.applied_ids = []
        self._applied_by_target.pop(getattr(result, "target", ""), None)
        self._refresh_after_edit(f"Removed {len(ids)} reaction(s) from the model.")

    def _save_designed_pathway(self) -> None:
        results = self.pathway_panel.all_results()
        if not results:
            return
        from .widgets.dialog_util import choose_save_path
        path = choose_save_path(self, "Save designed pathway", "designed_pathway.gsmpath",
                                "GSM pathway (*.gsmpath *.json)")
        if not path:
            return
        if not path.lower().endswith((".gsmpath", ".json")):
            path += ".gsmpath"
        import json
        data = {"pathways": [pathway_design.result_to_dict(r) for r in results]}
        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)
        except OSError as exc:
            QMessageBox.critical(self, "Could not save", str(exc))
            return
        self.status_label.setText(f"Saved designed pathway to {os.path.basename(path)}.")

    def _load_designed_pathway(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load designed pathway", "", "GSM pathway (*.gsmpath *.json)")
        if not path:
            return
        import json
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            results = [pathway_design.result_from_dict(d) for d in data.get("pathways", [])]
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Could not load", str(exc))
            return
        if not results:
            QMessageBox.information(self, "Load designed pathway", "No pathways in that file.")
            return
        self._pathway_results = results
        self.pathway_panel.show_results(results)
        self.status_label.setText(
            f"Loaded {len(results)} designed pathway(s) from {os.path.basename(path)}. "
            "Load the matching reaction database to add or draw them.")

    def _display_pathway(self) -> None:
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None:
            return
        ids = result.reaction_ids
        if not ids:
            QMessageBox.information(self, "Display reaction network",
                                    "This result has no heterologous reactions to display.")
            return

        def build():
            # Combine host + pathway (keeping database reaction ids so they can be
            # highlighted), then hand the model to the graph dialog.
            native_ids = {m.id for m in self.project.model.metabolites}
            combined = pathway_design.apply_pathway(
                self.project.model, db, ids, match_namespace=True, rename=False)
            # Anything in the combined model that the host did not already have came
            # from the database — apply_pathway reuses natives via cross-namespace
            # identity, so this difference is exactly the de-novo set to paint green.
            added = [m.id for m in combined.metabolites if m.id not in native_ids]
            return combined, added
        ok, built = run_busy(self, "Building pathway map…", build,
                             title="Display pathway", cancelable=True)
        if not ok:
            if not was_cancelled(built):
                QMessageBox.warning(self, "Could not display pathway", str(built))
            return
        combined, added = built
        from .dialogs.pathway_graph_dialog import PathwayGraphDialog
        present = [rid for rid in ids if combined.reactions.has_id(rid)]
        dlg = PathwayGraphDialog(self, combined, present, target=result.target,
                                 added_metabolites=added)
        dlg.exec()

    @staticmethod
    def _plain_compound_name(pretty: str) -> str:
        """'(E)-oct-2-enal (E_oct_2_enal_c)' -> '(E)-oct-2-enal' for an online lookup."""
        name = pretty.strip()
        if name.endswith(")") and " (" in name:
            name = name[:name.rfind(" (")]
        return name.strip()

    def _fetch_missing_chemistry(self, missing: list) -> None:
        """Fetch reactions around a compound no loaded database can produce.

        When a search fails because a precursor has no producer, that is a gap in the
        DATA, not the search — so the fix is to go and get the chemistry. Fetching
        around the MISSING compound (rather than the target) is what closes the chain:
        the target's own neighbourhood is already present, it is the branch feeding it
        that is absent.
        """
        if not missing:
            return
        names = [self._plain_compound_name(m) for m in missing]
        pick, ok = QInputDialog.getItem(
            self, "Fetch missing chemistry",
            "No loaded database can make these compounds. Fetch reactions "
            "around which one?", names, 0, False)
        if not ok or not pick.strip():
            return
        # Which chemistry source? KEGG and Rhea are open; MetaCyc needs a licence and
        # has no open bulk API, so it is offered only as a local-import note.
        source, ok = QInputDialog.getItem(
            self, "Fetch missing chemistry",
            f"Fetch reactions around “{pick}” from which source?",
            ["KEGG (open)", "Rhea (open, ChEBI-based)",
             "MetaCyc (requires a licence — how?)"], 0, False)
        if not ok:
            return
        if source.startswith("MetaCyc"):
            QMessageBox.information(
                self, "MetaCyc access",
                "MetaCyc/BioCyc reaction data is licence-restricted and has no open bulk "
                "API, so it cannot be fetched automatically.<br><br>"
                "If your institution has a BioCyc subscription, download the reactions "
                "flat-file for your organism from biocyc.org and load it via "
                "<b>Pathway Design ▸ Load database</b>. MetaCyc content also reaches you "
                "indirectly through the MetaNetX universal (Fetch online), which "
                "reconciles KEGG + MetaCyc + BiGG + Rhea + SEED.")
            return
        steps, ok = QInputDialog.getInt(
            self, "Fetch missing chemistry",
            f"Expansion steps around “{pick}”\n"
            "(0 = only reactions touching it; 1–2 = include its neighbours):",
            2, 0, 3)
        if not ok:
            return
        use_rhea = source.startswith("Rhea")
        src_name = "Rhea" if use_rhea else "KEGG"

        def work():
            builder = (databases.build_rhea_pathway_db if use_rhea
                       else databases.build_kegg_pathway_db)
            model, label, cache_path = builder(pick.strip(), expand_steps=steps)
            return model, label, cache_path

        self._run_fetch(
            f"Fetching {src_name} reactions around “{pick}”…", work,
            then=lambda: QMessageBox.information(
                self, "Fetched",
                f"{src_name} reactions around “{pick}” were added as a database.\n\n"
                "Tick it in the Reaction databases panel, click “Load selected "
                "databases”, then run Predict pathway again — the search will now "
                "have the chemistry it was missing."))

    def _diagnose_pathway_flux(self) -> None:
        """Explain WHY the selected route carries no flux, and what to do about it.

        Two very different causes look identical to the user (both show 0): the route
        may be genuinely blocked by a dead-end intermediate, or it may be perfectly
        capable and simply not rewarded by the current objective (FBA maximises
        biomass; making the product does not help the cell grow). FVA tells them apart.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None:
            return
        ids = result.reaction_ids
        if not ids:
            return

        def build():
            from ..core import pathway_diagnostics
            native = {m.id for m in self.project.model.metabolites}
            eng = pathway_design.apply_pathway(self.project.model, db, ids,
                                               rename=False, target=result.target)
            d = pathway_diagnostics.diagnose(eng, ids, target_id=result.target,
                                             native_ids=native)
            # Second-pass producibility: of the host precursors this route grounds on,
            # which cannot carry flux in the current medium? (Warning, not a filter.)
            try:
                from ..core import pathway_search
                prod = pathway_search.route_producibility_report(
                    self.project.model, db, ids)
            except Exception:  # noqa: BLE001
                prod = None
            return d, prod

        ok, built = run_busy(self, "Analysing pathway flux…", build,
                             title="Flux feasibility", cancelable=True)
        if not ok:
            if not was_cancelled(built):
                QMessageBox.warning(self, "Could not analyse", str(built))
            return
        diag, prod = built

        from ..core import pathway_diagnostics as pdg
        titles = {pdg.OK: "This pathway carries flux",
                  pdg.NOT_INCENTIVISED: "The pathway works — the objective is the issue",
                  pdg.BLOCKED: "The pathway is blocked",
                  pdg.INFEASIBLE: "The model could not be analysed"}
        parts = [f"<b>{diag.summary}</b>"]
        if diag.max_production == diag.max_production:
            parts.append(f"Maximum production if the product is the objective: "
                         f"<b>{diag.max_production:.4g}</b>.")
        if diag.verdict == pdg.BLOCKED:
            if diag.last_carrying:
                parts.append(f"Last step that can carry flux: <b>{diag.last_carrying}</b>.")
            parts.append(f"Bottleneck: <b>{diag.bottleneck}</b>.")
            if diag.blocking_metabolites:
                parts.append("Dead-end metabolite(s): <b>"
                             + ", ".join(diag.blocking_metabolites[:6]) + "</b>.")
        if diag.recommendation:
            parts.append("<br><u>What to do</u><br>" + diag.recommendation)
        # Second-pass precursor producibility (FVA warning, never a hard filter).
        if prod is not None:
            status = prod.get("status", "ok")
            nonprod = prod.get("nonproducible", [])
            if status in ("fva_failed", "no_evidence", "no_reactions"):
                parts.append("<br><u>Precursor producibility</u><br>"
                             "⚠ Could not verify whether the host can supply this route's "
                             "precursors (the flux-variability check did not complete under "
                             "the current conditions). Treat the flux result above as "
                             "unconfirmed on that point.")
            elif nonprod:
                names = ", ".join(f"<b>{nm}</b>" for _id, nm in nonprod[:6])
                parts.append("<br><u>Precursor producibility</u><br>"
                             "⚠ This route grounds on host compound(s) that carry <b>zero "
                             "production flux</b> under the current medium: " + names
                             + ". The route is still valid, but it will likely carry no "
                             "flux until you enable these precursors (change the medium, "
                             "or add a route to make them).")
        box = QMessageBox(self)
        box.setWindowTitle(titles.get(diag.verdict, "Pathway flux"))
        box.setIcon(QMessageBox.Information if diag.verdict in (pdg.OK, pdg.NOT_INCENTIVISED)
                    else QMessageBox.Warning)
        box.setTextFormat(Qt.RichText)
        box.setText("<br><br>".join(parts))
        box.exec()

    def _branching_analysis(self) -> None:
        """EA-MNE branching: which native reactions compete for this route's
        intermediates, how much yield that diversion costs, and offer the competitors as
        one-click knockdown candidates. A route is never a line through empty space."""
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        ids = result.reaction_ids

        def build():
            from ..core import pathway_diagnostics
            eng = pathway_design.apply_pathway(self.project.model, db, ids,
                                               rename=False, target=result.target)
            return pathway_diagnostics.analyse_branching(eng, ids, target_id=result.target)

        ok, br = run_busy(self, "Analysing competition for this route's intermediates…",
                          build, title="Branching & competition", cancelable=True)
        if not ok:
            if not was_cancelled(br):
                QMessageBox.warning(self, "Could not analyse", str(br))
            return
        parts = [f"<b>{br.summary}</b>"]
        ly, ny = br.linear_yield, br.network_yield
        if ly == ly and ny == ny:
            parts.append(f"Ideal (linear) yield: <b>{ly:.4g}</b> · "
                         f"realistic (network) yield: <b>{ny:.4g}</b>.")
            if br.yield_loss == br.yield_loss and br.yield_loss > 1e-6:
                parts.append(f"Branching costs about <b>{br.yield_loss * 100:.0f}%</b> "
                             "of the ideal yield.")
        if br.per_intermediate:
            rows = []
            for met, others in list(br.per_intermediate.items())[:8]:
                rows.append(f"&nbsp;&nbsp;<b>{met}</b> is also consumed by: "
                            + ", ".join(others[:6]) + (" …" if len(others) > 6 else ""))
            parts.append("<br><u>Competing reactions (knockdown candidates)</u><br>"
                         + "<br>".join(rows))
            parts.append("<br><i>Reducing or knocking out these reactions keeps more "
                         "carbon on the route to your product.</i>")
        elif ly == ly:
            parts.append("Nothing competes for this route's intermediates under the "
                         "current conditions — the linear yield is realistic.")
        box = QMessageBox(self)
        box.setWindowTitle("Branching & competition")
        box.setIcon(QMessageBox.Information)
        box.setTextFormat(Qt.RichText)
        box.setText("<br><br>".join(parts))
        box.exec()

    @staticmethod
    def _resolve_product_sink_impl(model, target_id):
        """The reaction whose maximisation measures export of ``target_id``.

        `apply_pathway` adds a transport to the extracellular compartment and puts the
        exchange on the *_e* twin, so looking only at the cytosolic metabolite finds no
        boundary and any capacity/strategy calculation silently returns zero. Follow the
        transport first, then fall back to a demand on the metabolite itself.
        """
        if not model.metabolites.has_id(target_id):
            return None
        met = model.metabolites.get_by_id(target_id)
        for r in met.reactions:
            if r.boundary:
                return r.id
        base = target_id.rsplit("_", 1)[0]
        for comp in ("e", "p"):
            twin = f"{base}_{comp}"
            if model.metabolites.has_id(twin):
                for r in model.metabolites.get_by_id(twin).reactions:
                    if r.boundary:
                        return r.id
        try:
            return model.add_boundary(met, type="demand").id
        except Exception:  # noqa: BLE001
            return None

    def _pathway_strategies(self) -> None:
        """FSEOF over/under-expression targets for the selected route (L9).

        Runs on the HOST-INTEGRATED model — the route applied to the host — so it works
        for rule-based (RetroRules) routes as well as database ones. Previously strategy
        analysis was only reachable after manually applying the pathway and switching to
        the Analysis tab, which meant the hardest targets never got a strategy scan.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        ids = result.reaction_ids

        def build():
            from ..core.analysis import strain_design as sd
            eng = pathway_design.apply_pathway(self.project.model, db, ids,
                                               rename=False, target=result.target)
            sink = self._resolve_product_sink_impl(eng, result.target)
            if sink is None:
                raise RuntimeError(
                    "Could not find or create an export reaction for "
                    f"'{result.target}', so there is nothing to enforce.")
            # Passing the route lets every hit be grouped as pathway / precursor supply /
            # competing sink rather than one flat list of hundreds of rows (VI.11).
            return sd.run_fseof(eng, sink, n_steps=8, pathway_reactions=ids), sink

        ok, built = run_busy(self, "Scanning for over/under-expression targets…",
                             build, title="Find strategies (FSEOF)", cancelable=True)
        if not ok:
            if not was_cancelled(built):
                QMessageBox.warning(self, "Could not scan", str(built))
            return
        res, sink = built
        table = res.table
        if table is None or not len(table):
            QMessageBox.information(
                self, "Find strategies (FSEOF)",
                "No clear over- or under-expression targets were found for this route.\n\n"
                + (res.note or ""))
            return
        # Hide the transport/exchange passengers by default; the user can ask for all.
        from ..core.analysis import strain_design as sd
        groups = sd.group_fseof_targets(res, hide_generic=True)
        shown = sum(len(g) for g in groups.values())
        hidden = len(table) - shown
        if hidden > 0:
            ans = QMessageBox.question(
                self, "Find strategies (FSEOF)",
                f"<b>{shown} engineering target(s)</b> found, grouped by role "
                f"(pathway, precursor supply, competing sink).<br><br>"
                f"{hidden} transport/exchange reaction(s) are hidden: they track overall "
                "turnover rather than being useful targets.<br><br>"
                "Show those as well?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ans == QMessageBox.Yes:
                groups = sd.group_fseof_targets(res, hide_generic=False)
        import pandas as pd
        parts = []
        for role, sub in groups.items():
            block = sub.copy()
            block.insert(0, "group", role)
            parts.append(block)
        out = pd.concat(parts, ignore_index=True) if parts else table
        self._show_table(out, f"FSEOF → {sink} · {res.note}")

    def _mdf_analysis(self) -> None:
        """Thermodynamic feasibility of the selected route via the Max-min Driving Force.

        A stoichiometrically valid route can still be thermodynamically impossible. The
        MDF optimises metabolite concentrations to make the least-favourable step as
        favourable as it can be: MDF > 0 means the route can run; MDF ≤ 0 names the step
        that is stuck no matter the concentrations.
        """
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None or not result.reaction_ids:
            return
        ids = result.reaction_ids
        from ..core import thermodynamics as thermo
        if not thermo.mdf_suite_enabled():
            QMessageBox.information(
                self, "Thermodynamics suite is disabled",
                "Thermodynamic (MDF) analysis is an optional suite and is currently "
                "switched off.<br><br>Enable it in <b>Settings ▸ Preferences ▸ Enable MDF "
                "Suite</b> — it needs a one-time ~1.34 GB data download.")
            return
        if not thermo.cache_present():
            QMessageBox.information(
                self, "Thermodynamics data missing",
                "The thermodynamics data has not been downloaded yet. Open "
                "<b>Settings ▸ Preferences</b> and download it, then try again.")
            return

        # The MDF depends entirely on the assumed pH, ionic strength and — most of all —
        # the concentration range it is allowed to optimise over. Show those settings
        # first and let the user change them; the analysis runs on the second click.
        from .dialogs.mdf_settings_dialog import MDFSettingsDialog
        settings_dlg = MDFSettingsDialog(
            self, metabolites=self._route_metabolites(db, ids),
            initial=getattr(self, "_mdf_settings", None))
        if settings_dlg.exec() != QDialog.Accepted:
            return
        self._mdf_settings = settings_dlg.settings()   # remembered for the next run
        kwargs = dict(self._mdf_settings)
        assumptions = settings_dlg.summary()

        def build():
            eng = pathway_design.apply_pathway(self.project.model, db, ids,
                                               rename=False, target=result.target)
            return thermo.analyse_pathway_mdf(eng, ids, **kwargs)

        # The first run also downloads a ~69 MB parameter file. Saying so turns an
        # apparent hang into an explained wait — and it only happens once per machine.
        busy_message = "Estimating ΔrG′ and the pathway driving force…"
        if not thermo.parameters_present():
            busy_message = (
                "Downloading the thermodynamic parameter file (~69 MB).\n"
                "This happens once; later analyses start immediately.")
        ok, res = run_busy(
            self, busy_message,
            build, title="Thermodynamics (MDF)", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                # An engine/data problem is not the same as a route with no ΔrG′°, and
                # the two need different advice. Point at the log either way: a packaged
                # app has no console, so that file is the only record of the traceback.
                from ..core import stdio_guard
                QMessageBox.warning(
                    self, "Thermodynamic analysis could not run",
                    f"{res}<br><br>"
                    "<span style='color:#5f6368'>This is a problem with the "
                    "thermodynamics engine or its data, not with your pathway. The full "
                    f"technical detail was written to<br><code>{stdio_guard.log_path()}"
                    "</code></span>")
            return
        if not res.dg0_prime:
            # Name the reactions and the compounds behind them: "no data" that does not
            # say WHICH data leaves the user with nowhere to go (VI.5).
            QMessageBox.warning(
                self, "Not enough data to analyse this route",
                "<b>None of this route's reactions could be scored</b>, so no driving "
                "force can be computed.<br><br>"
                + (res.data_warning() or "")
                + "<br><br>The identified reactions do not carry enough information in "
                "the loaded databases — typically a compound with no cross-reference and "
                "no structure, which is common for rule-generated chemistry.<br><br>"
                "<b>What this does and does not mean:</b> the route is <i>unscored</i>, "
                "not infeasible. Do not read the absence of a warning as a clean bill of "
                "health.")
            return
        if res.missing:
            QMessageBox.information(
                self, "Partial thermodynamic coverage", res.data_warning())
        verdict = ("<b>Thermodynamically feasible</b> at physiological concentrations."
                   if res.feasible else
                   "<b>Not feasible</b> as written — at least one step cannot be driven "
                   "forward at any concentration in range.")
        if res.single_reaction:
            parts = [verdict,
                     f"Driving force of this single reaction: <b>{res.mdf:.2f} kJ/mol</b>.",
                     "<span style='color:#8a6d00'>⚠ This route has only one reaction, so "
                     "there are no shared intermediates to balance: the optimiser simply "
                     "puts every substrate at its concentration ceiling and the product at "
                     "its floor. Read the value as “favourable”, <b>not</b> as a pathway "
                     "MDF comparable with multi-step routes.</span>"]
        else:
            parts = [verdict,
                     f"Max-min Driving Force: <b>{res.mdf:.2f} kJ/mol</b> "
                     "(positive is feasible; larger is more robust).",
                     f"Least-favourable step (bottleneck): <b>{res.bottleneck}</b> "
                     f"— ΔrG′ ≈ {res.dg_prime.get(res.bottleneck, float('nan')):.1f} kJ/mol."]
        rows = []
        for rid in ids:
            if rid in res.dg_prime:
                rows.append(f"&nbsp;&nbsp;<b>{rid}</b>: ΔrG′° "
                            f"{res.dg0_prime[rid]:+.1f} → ΔrG′ {res.dg_prime[rid]:+.1f} "
                            "kJ/mol")
        if rows:
            parts.append("<br><u>Per-reaction energies (standard → optimised)</u><br>"
                         + "<br>".join(rows))
        if res.missing:
            parts.append("<br><i>No ΔrG′° for: " + ", ".join(res.missing[:6])
                         + (" …" if len(res.missing) > 6 else "")
                         + " — these were left out of the driving-force calculation.</i>")
        # State the physiology the number rests on — the settings actually used, not a
        # fixed sentence: an MDF is only meaningful relative to those assumptions.
        parts.append(
            f"<br><span style='color:#5f6368'>Assumptions used: {assumptions} "
            "Feasibility calls close to zero are sensitive to these — re-run with a "
            "narrower concentration range to test how robust this verdict is.</span>")
        box = QMessageBox(self)
        box.setWindowTitle("Thermodynamics (MDF)")
        box.setIcon(QMessageBox.Information if res.feasible else QMessageBox.Warning)
        box.setTextFormat(Qt.RichText)
        box.setText("<br><br>".join(parts))
        box.exec()

    def _retry_for_flux_carrying_route(self) -> None:
        """Search again for an alternative route that can actually run, forbidding the
        step(s) that block the current one."""
        result = self.pathway_panel.current_result()
        if self.project is None or result is None or not result.reaction_ids:
            return
        forbid = set(result.blocked_by or result.reaction_ids)
        self._pathway_forbidden = getattr(self, "_pathway_forbidden", set()) | forbid
        QMessageBox.information(
            self, "Find a route that runs",
            "Searching again for an alternative route, avoiding "
            + ", ".join(sorted(forbid)[:4])
            + (" …" if len(forbid) > 4 else "")
            + ".\n\nIf no alternative exists, the database may simply not contain "
              "another way to make this compound.")
        self._predict_pathway(result.target, self.pathway_panel.min_flux.value(),
                              forbidden=self._pathway_forbidden)

    def _draw_pathway_scheme(self) -> None:
        result = self.pathway_panel.current_result()
        db = self._combined_included_model()
        if self.project is None or result is None or db is None:
            return
        ids = result.reaction_ids
        if not ids:
            QMessageBox.information(self, "Draw pathway scheme",
                                    "This result has no heterologous reactions to draw.")
            return
        from ..core.network_graph import short_metabolite_name
        from .widgets.structure_fetcher import fetch_structure_png, metabolite_structure_hints

        target = result.target

        def build():
            met_ids, steps = pathway_design.pathway_chain(db, ids, target)
            if not met_ids:      # couldn't linearise — fall back to reaction order
                met_ids, steps = pathway_design.pathway_chain(db, ids, ids[-1])
            nodes = []
            for mid in met_ids:
                met = db.metabolites.get_by_id(mid) if db.metabolites.has_id(mid) else None
                name = short_metabolite_name(mid, (met.name if met else "") or "")
                smiles = inchi = inchikey = kegg = chebi = ""
                if met is not None:
                    smiles, inchi, inchikey, kegg, chebi = metabolite_structure_hints(met)
                img = fetch_structure_png(name=(met.name if met else "") or name,
                                          inchikey=inchikey, smiles=smiles, inchi=inchi,
                                          kegg=kegg, chebi=chebi, size=160)
                nodes.append({"name": name, "img": img})
            arrows = []
            for rxn, consumed, produced, cosubs in steps:
                enzyme = rxn.name if (rxn.name and rxn.name != rxn.id) \
                    else pathway_design.readable_reaction_id(rxn)
                # Extra carbon-skeleton substrates (a reaction fusing two backbones)
                # are drawn as merge inputs with their own structure image.
                merges = []
                for c in cosubs:
                    m = db.metabolites.get_by_id(c) if db.metabolites.has_id(c) else None
                    mname = short_metabolite_name(c, (m.name if m else "") or "")
                    ms = mi = mk = mkegg = mchebi = ""
                    if m is not None:
                        ms, mi, mk, mkegg, mchebi = metabolite_structure_hints(m)
                    merges.append({"name": mname,
                                   "img": fetch_structure_png(
                                       name=(m.name if m else "") or mname,
                                       inchikey=mk, smiles=ms, inchi=mi,
                                       kegg=mkegg, chebi=mchebi, size=90)})
                arrows.append({
                    "enzyme": enzyme,
                    "consumed": ", ".join(short_metabolite_name(c, "") for c in consumed),
                    "produced": ", ".join(short_metabolite_name(c, "") for c in produced),
                    "merges": merges,
                })
            return nodes, arrows

        ok, res = run_busy(self, "Drawing pathway scheme (fetching structures)…", build,
                           title="Pathway scheme", cancelable=True)
        if not ok:
            if not was_cancelled(res):
                QMessageBox.warning(self, "Could not draw pathway scheme", str(res))
            return
        nodes, arrows = res
        if not nodes:
            QMessageBox.information(self, "Draw pathway scheme",
                                    "Could not lay out this pathway as a linear scheme.")
            return
        from .dialogs.pathway_scheme_dialog import PathwaySchemeDialog
        PathwaySchemeDialog(self, nodes, arrows, target=target).exec()

    def _explore_alternative_pathway(self) -> None:
        db = self._combined_included_model()
        if self.project is None or db is None:
            return
        target = self._resolve_target_in(db, self.pathway_panel.current_target())
        if not target:
            return
        forbidden = self.pathway_panel.all_reaction_ids()
        host = self.project.model.copy()

        def finish(results):
            new = results[0] if results else None
            if new is None or new.reactions.empty:
                QMessageBox.information(
                    self, "Explore alternative pathways",
                    "No further alternative pathway was found. The routes already shown may be "
                    "the only ones the loaded databases support — try increasing the pathway "
                    "length or including more databases.")
                return
            self.pathway_panel.append_result(new)
            self.status_label.setText(f"Found an alternative pathway to {target}.")

        ok, res = run_busy(
            self, f"Exploring an alternative pathway to {target}…",
            lambda: pathway_design.find_pathways(
                host, target, db, n_alternatives=1, forbidden_reactions=set(forbidden)),
            title="Pathway design", after=finish, after_message="Preparing results…",
            cancelable=True)
        if not ok and not was_cancelled(res):
            QMessageBox.critical(self, "Explore alternative pathways failed", str(res))

    def _run_efm(self, cfg: dict) -> None:
        cat_name = self.analysis_panel.current_category()
        if not cat_name or not self.project.categories.has(cat_name):
            QMessageBox.information(
                self, "Elementary flux modes",
                "EFM analysis runs on a small subnetwork. Choose a category in the "
                "Analysis 'Scope' selector first (create one in the Categories panel).")
            return
        ids = set(self.project.categories.get(cat_name).reaction_ids)
        if not ids:
            QMessageBox.information(self, "Elementary flux modes", "The category is empty.")
            return
        model = self.project.model.copy()
        max_modes = cfg.get("max_modes", 200)

        def done(result):
            if result.table.empty:
                QMessageBox.information(self, "Elementary flux modes", result.note)
                return
            self._show_table(result.table, f"Elementary flux modes · {cat_name} · {result.note}")
        self._run_job(lambda: pathways.enumerate_efms(model, ids, max_modes=max_modes), done,
                      title=f"Elementary flux modes · {cat_name}", kind="efm")

    def _edit_consortia_objective(self) -> bool:
        """Open the consortia objective settings (dominance + no-starvation floor) and
        store them in the project. Returns True if the user applied changes."""
        members = self.project.settings.get("community_members") if self.project else None
        if not members:
            QMessageBox.information(
                self, "Consortia objective",
                "This is not a community model. Build one from Tools ▸ Consortia "
                "modelling ▸ Build Community Model first.")
            return False
        from .dialogs.consortia_objective_dialog import ConsortiaObjectiveDialog
        prev = self.project.settings.get("consortia_objective", {})
        dlg = ConsortiaObjectiveDialog(
            self, list(members.keys()), weights=prev.get("weights"),
            min_growth_fraction=prev.get("min_growth_fraction", 0.0))
        if dlg.exec() != ConsortiaObjectiveDialog.Accepted:
            return False
        self.project.settings["consortia_objective"] = dlg.values()
        return True

    def _run_community_growth(self, cfg: dict) -> None:
        members = self.project.settings.get("community_members")
        if not members:
            QMessageBox.information(
                self, "Community member growth",
                "No community model loaded. Use Tools ▸ Consortia modelling ▸ "
                "Build Community Model first.")
            return
        # Let the user set dominance + the no-starvation floor before running (#F2).
        if not self._edit_consortia_objective():
            return
        opts = self.project.settings.get("consortia_objective", {})
        weights = opts.get("weights")
        min_growth = float(opts.get("min_growth_fraction", 0.0) or 0.0)
        model = self.project.model.copy()

        def done(df):
            if df.attrs.get("infeasible") or (
                    "growth" in df.columns and df["growth"][:-1].isna().all()):
                QMessageBox.warning(
                    self, "Community member growth",
                    df.attrs.get("infeasible")
                    or "No feasible community state at these settings — lower the "
                       "minimum-growth requirement and re-run.")
            self._show_table(df, "Community member growth (FBA)")
        self._run_job(
            lambda: community_core.member_growth_table(
                model, members, weights=weights, min_growth_fraction=min_growth),
            done, title="Community member growth (FBA)", kind="generic")

    # ----- Phase 5: omics & energy -------------------------------------
    def _ensure_omics_tab(self):
        """Create the Omics tab lazily and make it visible (#F4). The tab exists
        only once an omics dataset has been prepared/loaded — a mapping to the model."""
        if self.omics_panel is None:
            from .panels.omics_panel import OmicsPanel
            self.omics_panel = OmicsPanel()
        if self.tabs.indexOf(self.omics_panel) < 0:
            idx = self.tabs.indexOf(self.growth_panel)
            self.tabs.insertTab(idx + 1 if idx >= 0 else self.tabs.count(),
                                self.omics_panel, "Omics")
        return self.omics_panel

    def _hide_omics_tab(self) -> None:
        if self.omics_panel is not None:
            idx = self.tabs.indexOf(self.omics_panel)
            if idx >= 0:
                self.tabs.removeTab(idx)
            self.omics_panel.clear_all()

    def _refresh_omics_tab(self) -> None:
        """Show/populate or hide the Omics tab based on the current project's datasets."""
        prepared = {} if self.project is None else self.project.datasets.get("omics_prepared", {})
        if not prepared:
            self._hide_omics_tab()
            return
        panel = self._ensure_omics_tab()
        panel.clear_all()
        for name, d in prepared.items():
            panel.add_dataset(name, d.get("values", {}), summary=d.get("summary", ""),
                              kind=d.get("kind", ""))

    def _prepare_omics_dataset(self) -> None:
        from .dialogs.omics_dataset_dialog import OmicsDatasetDialog
        dlg = OmicsDatasetDialog(self, self.project.model)
        if dlg.exec() != OmicsDatasetDialog.Accepted or dlg.dataset is None:
            return
        ds = dlg.dataset
        from ..core import omics_prep as op
        name = f"{ds.kind} ({ds.summary.n_model_targets} mapped)"
        store = self.project.datasets.setdefault("omics_prepared", {})
        store[name] = {"values": ds.values, "summary": ds.summary.text(), "kind": ds.kind}
        # Gene-based data feeds eFlux/GIMME directly; metabolomics is kept for reference.
        if ds.kind in op.GENE_KINDS:
            self.project.datasets["expression"] = ds.values
        else:
            self.project.datasets["metabolomics"] = ds.values
        self.project._dirty = True
        panel = self._ensure_omics_tab()
        panel.add_dataset(name, ds.values, summary=ds.summary.text(), kind=ds.kind)
        self.tabs.setCurrentWidget(panel)
        follow = ("You can now run eFlux or GIMME from the Omics & energy section."
                  if ds.kind in op.GENE_KINDS else
                  "Metabolomics values are mapped to model metabolites for reference.")
        QMessageBox.information(
            self, "Omics dataset prepared",
            f"Mapped {ds.summary.n_model_targets} of {ds.summary.n_model_total} "
            f"model {'genes' if ds.kind in op.GENE_KINDS else 'metabolites'} "
            f"({ds.summary.coverage*100:.0f}% coverage). {follow}")

    def _load_expression(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open gene-expression table", "",
            "Expression table (*.csv *.tsv *.txt);;All files (*)")
        if not path:
            return
        try:
            expr = omics.load_expression(path)
        except omics.OmicsError as exc:
            QMessageBox.critical(self, "Could not load expression data", str(exc))
            return
        self.project.datasets["expression"] = expr
        # How many of these ids actually match model genes (a quick coverage hint).
        model_genes = {g.id for g in self.project.model.genes}
        matched = sum(1 for k in expr if str(k) in model_genes)
        name = f"expression ({matched} matched)"
        summary = (f"Loaded a ready expression table: {len(expr)} values, "
                   f"{matched} match model genes ({matched/max(len(model_genes),1)*100:.0f}% "
                   "coverage). If this is low, use 'Prepare omics dataset…' to map identifiers.")
        self.project.datasets.setdefault("omics_prepared", {})[name] = {
            "values": expr, "summary": summary, "kind": "transcriptomics"}
        self.project._dirty = True
        panel = self._ensure_omics_tab()
        panel.add_dataset(name, expr, summary=summary, kind="transcriptomics")
        QMessageBox.information(
            self, "Expression data loaded",
            f"Loaded expression values for {len(expr)} genes ({matched} match model genes). "
            "You can now run eFlux or GIMME from the Omics & energy section.")
        self.status_label.setText(f"Loaded expression data ({len(expr)} genes).")

    def _expression_or_warn(self):
        expr = self.project.datasets.get("expression")
        if not expr:
            QMessageBox.information(
                self, "No expression data",
                "Load a gene-expression table first (Omics & energy ▸ Load expression data).")
            return None
        return expr

    def _run_eflux(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        expr = self.project.datasets.get("expression")

        def done(result):
            self._show_table(result.table, f"eFlux · objective = {result.objective_value:.6g} · {scope}")
        self._run_job(lambda: omics.run_eflux(model, expr), done,
                      title=f"eFlux · {scope}", kind="eflux")

    def _run_gimme(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        expr = self.project.datasets.get("expression")

        def done(result):
            self._show_table(result.table, f"GIMME · {result.note}")
        self._run_job(
            lambda: omics.run_gimme(model, expr, threshold=cfg["threshold"],
                                    objective_fraction=cfg["objective_fraction"]), done,
            title=f"GIMME · {scope}", kind="gimme")

    def _run_atpm(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        atpm_id = cfg["atpm_id"].strip()
        if not atpm_id:
            return

        def done(result):
            self._show_table(result.table, f"ATP maintenance sensitivity · {scope}")
            self.analysis_panel.plot_view.plot_line(
                result.table, "atp_maintenance", "growth",
                title="ATP maintenance sensitivity", xlabel="ATP maintenance flux",
                ylabel="Growth rate")
            self.analysis_panel.result_tabs.setCurrentWidget(self.analysis_panel.plot_view)
        self._run_job(
            lambda: omics.run_atpm_sensitivity(model, atpm_id=atpm_id, points=cfg["points"]), done,
            title=f"ATP maintenance sensitivity · {scope}", kind="atpm_sensitivity")

    # ----- Phase 4: strain design --------------------------------------
    _KO_METHODS = {
        "OptKnock": "optknock", "RobustKnock": "robustknock",
        "OptCouple": "optcouple", "Heuristic (evolutionary)": "heuristic",
    }
    _KO_APPROACHES = {"Best (optimal)": "best", "Any (faster)": "any", "Diverse set": "populate"}

    def _warn_if_not_exchange(self, model, product: str, title: str) -> bool:
        """C15: nudge the user toward an exchange reaction as the maximization target.
        Returns True to proceed, False to cancel."""
        rxn = model.reactions.get_by_id(product) if model.reactions.has_id(product) else None
        if rxn is not None and not rxn.boundary:
            return QMessageBox.question(
                self, title,
                f"'{product}' is an internal reaction, not an exchange.\n\n"
                "For product design you usually want to maximize the EXCHANGE reaction that "
                "secretes your product (e.g. EX_<product>_e) so growth is coupled to real "
                "export. Maximizing an internal step can be misleading.\n\nContinue anyway?",
                QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes
        return True

    def _run_knockout_design(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        product = cfg["product"]
        if not model.reactions.has_id(product):
            QMessageBox.warning(self, "Knockout strain design", f"No reaction '{product}'.")
            return
        if not self._warn_if_not_exchange(model, product, "Knockout strain design"):
            return
        method = self._KO_METHODS.get(cfg["method"], "optknock")
        solver = "auto" if cfg["solver"] == "Auto" else cfg["solver"]
        approach = self._KO_APPROACHES.get(cfg["approach"], "best")
        self.status_label.setText(f"Running {cfg['method']} — please wait…")

        def done(result):
            if result.table.empty:
                QMessageBox.information(self, "Knockout strain design", result.note)
                return
            self._show_table(result.table, f"{result.method} → {product} · {result.note}")
        self._run_job(
            lambda: sd.run_knockout_design(
                model, product, method=method, solver=solver, approach=approach,
                biomass_id=cfg["biomass"] or None,
                max_knockouts=cfg["max_knockouts"], max_solutions=cfg["max_solutions"],
                min_growth_fraction=cfg["min_growth_fraction"],
                exclude_exchanges=cfg["exclude_exchanges"],
                exclude_transport=cfg["exclude_transport"],
                exclude_blocked=cfg["exclude_blocked"],
                time_limit=cfg["time_limit"],
                population=cfg["population"], generations=cfg["generations"]),
            done, title=f"{cfg['method']} strain design → {product}", kind="knockout")

    def _run_fseof(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        product = cfg["product"]
        if not model.reactions.has_id(product):
            QMessageBox.warning(self, "FSEOF", f"No reaction '{product}'.")
            return
        if not self._warn_if_not_exchange(model, product, "FSEOF"):
            return

        def done(result):
            if result.table.empty:
                QMessageBox.information(self, "FSEOF", "No clear over/under-expression targets found.")
                return
            self._show_table(result.table, f"FSEOF → {product} · {result.note}")
        self._run_job(
            lambda: sd.run_fseof(model, product, n_steps=cfg["n_steps"],
                                 tolerance=cfg.get("tolerance", 1e-6),
                                 include_trending=cfg.get("include_trending", True)),
            done, title=f"FSEOF → {product}", kind="fseof")

    def _run_overproduction(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        target = (cfg.get("metabolite") or "").strip()
        if not target or not model.metabolites.has_id(target):
            QMessageBox.warning(self, "Overproduction", f"Choose a valid metabolite (got '{target}').")
            return

        def done(result):
            if result.table.empty:
                QMessageBox.information(
                    self, "Overproduction",
                    "No clear amplification / knock-down targets were found for "
                    f"{target}. The model may already be at its production limit.")
                return
            self._show_table(result.table, f"Overproduce {target} · {result.note}")
        self._run_job(
            lambda: sd.run_metabolite_overproduction(
                model, target, n_steps=cfg["n_steps"],
                tolerance=cfg.get("tolerance", 1e-6),
                include_trending=cfg.get("include_trending", True),
                secretion=cfg.get("secretion", False),
                gas=cfg.get("gas", False)), done,
            title=f"Overproduction targets for {target}", kind="overproduction")

    # ----- Phase 3: QC & gap-filling -----------------------------------
    def _run_quality_report(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        include_blocked = cfg.get("include_blocked", True)

        def done(report):
            self.analysis_panel.show_table(report.summary, f"Model quality report · {scope}")
        self._run_job(lambda: qc.quality_report(model, include_blocked=include_blocked), done,
                      title=f"Model quality report · {scope}", kind="quality_report")

    def _run_blocked_reactions(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        def done(report):
            df = report.blocked if not report.blocked.empty else __import__("pandas").DataFrame(
                {"reaction": ["(none — all reactions can carry flux)"]})
            self.analysis_panel.show_table(df, f"Blocked reactions · {scope}")
        self._run_job(lambda: qc.quality_report(model, include_blocked=True), done,
                      title=f"Finding blocked reactions · {scope}", kind="blocked_reactions")

    def _load_universal_model(self):
        """Prompt for a universal model file used as the gap-filling reaction pool."""
        QMessageBox.information(
            self, "Universal model needed",
            "Gap-filling needs a 'universal' model — a pool of candidate reactions to "
            "draw from (e.g. a BiGG universal model). Please choose that SBML/JSON file next.")
        path, _ = QFileDialog.getOpenFileName(self, "Open universal model", "", _MODEL_FILTER)
        if not path:
            return None
        try:
            return io_models.load_model(path)
        except io_models.ModelLoadError as exc:
            QMessageBox.critical(self, "Could not open universal model", str(exc))
            return None

    def _run_gapfill_growth(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        universal = self._load_universal_model()
        if universal is None:
            return
        import pandas as pd
        lb = cfg.get("lower_bound", 0.05)

        def done(solutions):
            rows = [{"solution": i + 1, "reactions_to_add": ", ".join(s)}
                    for i, s in enumerate(solutions)] or [{"solution": "—", "reactions_to_add": "no solution found"}]
            self.analysis_panel.show_table(pd.DataFrame(rows), f"Gap-fill for growth · {scope}")
        self._run_job(lambda: gapfill.gapfill_for_growth(model, universal, lower_bound=lb), done,
                      title=f"Gap-filling for growth · {scope}", kind="gapfill_growth")

    def _run_gapfill_metabolite(self, analysis_id: str, model, scope: str, cfg: dict) -> None:
        target = cfg.get("metabolite", "").strip()
        if not target or not model.metabolites.has_id(target):
            QMessageBox.warning(self, "Gap-fill", f"Choose a valid metabolite (got '{target}').")
            return
        universal = self._load_universal_model()
        if universal is None:
            return
        import pandas as pd
        lb = cfg.get("lower_bound", 0.05)

        def done(solutions):
            rows = [{"solution": i + 1, "reactions_to_add": ", ".join(s)}
                    for i, s in enumerate(solutions)] or [{"solution": "—", "reactions_to_add": "no solution found"}]
            self.analysis_panel.show_table(pd.DataFrame(rows), f"Gap-fill to produce {target} · {scope}")
        self._run_job(lambda: gapfill.gapfill_for_metabolite(model, universal, target, lower_bound=lb), done,
                      title=f"Gap-filling for {target} · {scope}", kind="gapfill_metabolite")

    def _pick_metabolite(self, model, title: str, prompt: str) -> Optional[str]:
        ids = [m.id for m in model.metabolites]
        if not ids:
            return None
        choice, ok = QInputDialog.getItem(self, title, prompt, ids, 0, True)
        return choice if ok and choice else None

    def _pick_reaction(self, model, title: str, prompt: str) -> Optional[str]:
        """Let the user choose a reaction id, defaulting to the explorer selection."""
        ids = [r.id for r in model.reactions]
        selected = self.explorer.selected_reaction()
        default = 0
        if selected is not None and selected.id in ids:
            default = ids.index(selected.id)
        choice, ok = QInputDialog.getItem(self, title, prompt, ids, default, True)
        return choice if ok and choice else None

    # ----- job plumbing -------------------------------------------------
    def _run_job(self, fn, on_done, *, title: str = "Working…", kind: str = "generic") -> None:
        """Run ``fn`` as a tracked background job. ``title`` is the specific
        description shown in the status bar; ``kind`` keys the time estimate.
        Multiple jobs may run at once — the UI stays usable throughout."""
        def done(result):
            try:
                on_done(result)
            except Exception as exc:  # noqa: BLE001 - never let rendering crash the app
                QMessageBox.critical(self, "Could not display results", str(exc))
        self.jobs.submit(fn, title=title, kind=kind, on_done=done, on_error=self._job_failed)

    def _job_failed(self, message: str) -> None:
        QMessageBox.critical(self, "Analysis failed", message)

    def _on_fba_done(self, method, result, whole_model, want_shadow, scope) -> None:
        if not result.is_optimal:
            QMessageBox.warning(
                self, f"{method} not optimal",
                f"The solver returned status '{result.status}'. The model may be "
                "over-constrained — check the growth medium and reaction bounds.")
            self.objective_label.setText(f"{method}: {result.status}")
            return
        self._last_result = result
        if want_shadow:
            sp = result.shadow_prices
            if sp is None or sp.empty:
                QMessageBox.information(self, "Shadow prices",
                                        "No shadow prices were returned for this solution.")
                return
            df = sp.rename("shadow_price").to_frame()
            df.index.name = "metabolite"
            self._show_table(df.reset_index(), f"Shadow prices · objective {result.objective_value:.6g}")
        else:
            note = getattr(result, "note", "")
            header = (f"{method} · objective = {result.objective_value:.6g} · "
                      f"{len(result.active_fluxes())} reactions carry flux · {scope}")
            if note:
                header += f" · ⚠ {note}"
            self._show_table(result.flux_table(), header)
            if note:
                QMessageBox.information(self, f"{method} note", note)
            if whole_model:
                fluxes = result.fluxes.to_dict()
                self.explorer.set_fluxes(fluxes)
                self.network_view.set_fluxes(fluxes)
        self.info.update_model(self.project.model, self.project.diff(), result.objective_value)
        self.objective_label.setText(f"{method} objective: {result.objective_value:.6g}")
        self.status_label.setText(f"{method} complete ({scope}).")

    def _open_manage_data(self) -> None:
        """Open the cache manager (databases + molecule images) as a pop-up window."""
        dlg = QDialog(self)
        dlg.setWindowTitle("Manage Data — cached databases & images")
        dlg.resize(720, 560)
        lay = QVBoxLayout(dlg)
        panel = SettingsPanel()
        panel.refresh()
        lay.addWidget(panel)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        dlg.exec()

    def _set_results_dir(self) -> None:
        from PySide6.QtWidgets import QFileDialog
        start = self._results_dir or os.path.expanduser("~")
        folder = QFileDialog.getExistingDirectory(self, "Results output folder", start)
        if folder:
            self._results_dir = folder
            self.status_label.setText(f"Analysis results will be saved to: {folder}")

    def _ensure_results_dir(self) -> Optional[str]:
        """Return the results folder, prompting the user to choose one if unset."""
        if not self._results_dir:
            QMessageBox.information(
                self, "Results folder",
                "No results folder has been set yet. Choose a folder to save results into.")
            self._set_results_dir()
        return self._results_dir

    def _result_filename(self, base: str, ext: str) -> str:
        import re
        import time as _t
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", (base or "result").split("·")[0].strip())[:60]
        slug = slug or "result"
        return os.path.join(self._results_dir,
                            f"{slug}_{_t.strftime('%Y%m%d_%H%M%S')}.{ext}")

    def _save_analysis_results(self) -> None:
        kind, payload, title = self.analysis_panel.active_result()
        if kind is None:
            QMessageBox.information(self, "Save Results",
                                    "There are no results to save yet. Run an analysis first.")
            return
        if not self._ensure_results_dir():
            return
        saved: list = []
        try:
            if kind == "table":
                path = self._result_filename(title, "csv")
                payload.to_csv(path, index=False)
                saved.append(path)
            else:  # plot
                png = self._result_filename(title, "png")
                pdf = self._result_filename(title, "pdf")
                payload.save(png)
                payload.save(pdf)
                saved.extend([png, pdf])
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Could not save results", str(exc))
            return
        names = ", ".join(os.path.basename(p) for p in saved)
        self.status_label.setText(f"Saved results to {self._results_dir}: {names}")
        QMessageBox.information(self, "Results saved",
                                f"Saved to your results folder:\n\n{chr(10).join(saved)}")

    def _save_pathway_results(self) -> None:
        result = self.pathway_panel.current_result()
        if result is None or result.reactions.empty:
            QMessageBox.information(self, "Save Results", "There is no pathway to save yet.")
            return
        if not self._ensure_results_dir():
            return
        try:
            path = self._result_filename(f"pathway_{result.target}", "csv")
            result.reactions.to_csv(path, index=False)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Could not save results", str(exc))
            return
        self.status_label.setText(f"Saved pathway to {path}")
        QMessageBox.information(self, "Results saved", f"Saved the pathway to:\n\n{path}")

    def _autosave_result(self, df, header: str) -> None:
        if not self._results_dir or df is None or getattr(df, "empty", True):
            return
        import re
        import time as _t
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", header.split("·")[0].strip())[:60] or "result"
        path = os.path.join(self._results_dir, f"{slug}_{_t.strftime('%Y%m%d_%H%M%S')}.csv")
        try:
            df.to_csv(path, index=False)
        except Exception:  # noqa: BLE001 - saving is a convenience, never fatal
            pass

    def _show_table(self, df, header: str) -> None:
        df = self._filter_excluded(df)
        self.analysis_panel.show_table(df, header)
        self._last_table = df          # kept for the Plot Gallery
        self.tabs.setCurrentWidget(self.analysis_panel)
        self.status_label.setText(header)
        self._autosave_result(df, header)
        # Offer a flux-network view when the table carries per-reaction fluxes.
        self._flux_data = self._extract_flux_data(df)
        self.analysis_panel.set_flux_available(self._flux_data is not None)

    @staticmethod
    def _extract_flux_data(df):
        """From a result table, build ``{reaction_id: (flux, label)}`` when it holds
        per-reaction fluxes (FBA/pFBA/FVA/FSEOF/MOMA/ROOM…). Returns None otherwise."""
        if df is None or getattr(df, "empty", True):
            return None
        cols = {c.lower(): c for c in df.columns}
        rcol = next((cols[c] for c in ("reaction", "reaction_id", "id", "rxn") if c in cols), None)
        if rcol is None:
            return None
        fcol = next((cols[c] for c in ("flux", "fluxes", "value", "flux_wt") if c in cols), None)
        lo = next((cols[c] for c in ("minimum", "flux_minimum", "lower") if c in cols), None)
        hi = next((cols[c] for c in ("maximum", "flux_maximum", "upper") if c in cols), None)
        if fcol is None and not (lo and hi):
            return None
        out = {}
        for _, row in df.iterrows():
            rid = str(row[rcol])
            try:
                if fcol is not None:
                    val = float(row[fcol])
                    label = f"{val:.3g}"
                    if lo and hi:
                        label += f" [{float(row[lo]):.2g}, {float(row[hi]):.2g}]"
                else:
                    lov, hiv = float(row[lo]), float(row[hi])
                    val = (lov + hiv) / 2.0
                    label = f"[{lov:.2g}, {hiv:.2g}]"
            except (TypeError, ValueError):
                continue
            out[rid] = (val, label)
        return out or None

    def _display_fluxes(self) -> None:
        if self.project is None or not getattr(self, "_flux_data", None):
            QMessageBox.information(self, "Display fluxes",
                                    "Run a flux analysis (FBA, pFBA, FVA, FSEOF…) first.")
            return
        model = self.project.model
        # Scope dialog: whole model / a compartment / a category.
        dlg = QDialog(self)
        dlg.setWindowTitle("Display fluxes in network")
        form = QVBoxLayout(dlg)
        form.addWidget(QLabel("Which reactions to draw (only non-zero fluxes are shown):"))
        combo = QComboBox()
        combo.addItem("All reactions carrying flux", ("all", None))
        for comp in sorted(getattr(model, "compartments", {}) or {}):
            combo.addItem(f"Compartment: {comp}", ("comp", comp))
        for name in self.project.categories.names():
            combo.addItem(f"Category: {name}", ("cat", name))
        form.addWidget(combo)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, dlg)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        form.addWidget(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        kind, arg = combo.currentData()

        flux = {r: v for r, (v, _t) in self._flux_data.items()}
        labels = {r: t for r, (_v, t) in self._flux_data.items()}
        # Reactions to include: those with (near-)non-zero flux within the chosen scope.
        active = {r for r, v in flux.items()
                  if abs(v) > 1e-9 and model.reactions.has_id(r)}
        if kind == "comp":
            active = {r for r in active
                      if arg in (model.reactions.get_by_id(r).compartments or set())}
        elif kind == "cat" and self.project.categories.has(arg):
            catset = set(self.project.categories.get(arg).reaction_ids)
            active = {r for r in active if r in catset}
        if not active:
            QMessageBox.information(self, "Display fluxes",
                                    "No reactions carry flux in the selected scope.")
            return
        from .dialogs.flux_map_dialog import FluxMapDialog
        FluxMapDialog(self, model, sorted(active),
                      {r: flux[r] for r in active}, {r: labels[r] for r in active}).exec()

    def _filter_excluded(self, df):
        """Drop transport/exchange reaction rows from a result table when the
        analysis's settings asked to exclude them."""
        excluded = getattr(self, "_excluded_rxn_ids", None)
        if not excluded or df is None or getattr(df, "empty", True):
            return df
        if "reaction" in df.columns:
            return df[~df["reaction"].isin(excluded)].reset_index(drop=True)
        return df

    # ----- selection relays --------------------------------------------
    def _on_reaction_selected(self, rxn) -> None:
        flux = None
        if self._last_result is not None and self._last_result.is_optimal:
            flux = self._last_result.fluxes.get(rxn.id)
        self.info.show_reaction(rxn, flux)

    # ----- refresh helpers ---------------------------------------------
    def _refresh_after_edit(self, message: str) -> None:
        self._last_result = None
        self._refresh_all()
        self.status_label.setText(message)

    def _update_action_availability(self) -> None:
        """Grey out everything that needs a model when none is loaded.

        This is the real fix behind "menu entries that do nothing on click". Handlers
        like ``edit_medium`` open with ``if self.project is None: return`` — fully
        implemented, but silently inert without a model, so the menu accepted clicks and
        ignored them. Disabling the actions makes the interface state visible instead of
        leaving the user to guess whether a feature is broken or simply unavailable.
        """
        has_model = self.project is not None
        # Reactions, medium, save/export and the undo pair are already handled — undo
        # and redo also depend on the history, which only that method knows about.
        self._update_action_states()
        for attr in ("act_model_info", "act_community"):
            action = getattr(self, attr, None)
            if action is not None:
                action.setEnabled(has_model)
        for menu in getattr(self, "_tools_model_menus", ()):
            menu.setEnabled(has_model)
        button = getattr(self, "_network_button", None)
        if button is not None:
            button.setEnabled(has_model)

    def _refresh_all(self) -> None:
        # Repopulating every panel resizes a lot of widgets at once; without this the
        # window drops out of Maximized as a side effect of loading a model.
        with self.hold_window_state():
            self._refresh_all_inner()

    def _refresh_all_inner(self) -> None:
        self._update_action_availability()
        self._sync_regulation_status()
        if self.project is None:
            return
        model = self.project.model
        self.explorer.set_model(model)
        # The dynamic panel builds its condition list and substrate table from this
        # model's own exchanges, so it has to be told when the model changes.
        self.dynamic_panel.set_model(model)
        self.network_view.set_added_metabolites(self.project.diff().get("added_metabolites", []))
        self.network_view.set_model(model)
        try:
            from ..core import physiology
            self.analysis_panel.set_physiology(physiology.summarize(model))
        except Exception:  # noqa: BLE001 - a banner must never break a refresh
            self.analysis_panel.set_physiology(None)
        self.growth_panel.set_model(model)
        self.strategy_explorer.set_model(model)
        self.strategy_explorer.set_strategies(self.project.strategies)
        self.escher_explorer.set_model(model)
        self.escher_explorer.set_strategies(self.project.strategies)
        self.info.update_model(model, self.project.diff())
        self.info.clear_selection()
        self.analysis_panel.objective_bar.set_model(model)
        self.analysis_panel.objective_bar.set_community(
            bool(self.project.settings.get("community_members")))
        self._refresh_omics_tab()
        self.pathway_panel.set_model(model)
        self._refresh_pathway_targets()
        self._refresh_categories()
        self.objective_label.setText("")
        self._update_title()
        self._update_action_states()

    def _update_title(self) -> None:
        title = __app_name__
        if self.project is not None:
            name = (self.project.project_path or self.project.source_path or "untitled")
            mark = "*" if self.project.is_modified else ""
            title = f"{__app_name__} — {os.path.basename(name)}{mark}"
        self.setWindowTitle(title)

    def _update_action_states(self) -> None:
        has = self.project is not None
        for a in (self.act_save_project, self.act_save_project_as, self.act_export_model,
                  self.act_add_rxn, self.act_edit_rxn, self.act_remove_rxn, self.act_set_obj,
                  self.act_media):
            a.setEnabled(has)
        self.act_undo.setEnabled(has and self.project.can_undo)
        self.act_redo.setEnabled(has and self.project.can_redo)

    def _on_result_context(self, cells: list, global_pos) -> None:
        """Right-click on a results-table row: actions for the reaction it names."""
        if self.project is None:
            return
        model = self.project.model
        rid = next((c for c in cells if model.reactions.has_id(c)), None)
        if rid is None:
            return
        rxn = model.reactions.get_by_id(rid)
        menu = QMenu(self)
        show_menu = menu.addMenu("Show in map")
        for steps in (1, 2, 3):
            act = show_menu.addAction(f"{steps} step{'s' if steps > 1 else ''} around")
            act.triggered.connect(lambda _=False, s=steps, r=rid: self._show_in_map(r, s))
        menu.addAction("Reaction information…", lambda: self._show_reaction_details(rxn))
        menu.addAction("Edit reaction…", lambda: self._edit_reaction(rxn))
        menu.addAction("Set as objective", lambda: self._set_objective_id(rid))
        menu.addSeparator()
        menu.addAction("Add to category…", lambda: self._add_reactions_to_category([rid]))
        menu.exec(global_pos)

    def _suggest_enzymes(self, rxn) -> None:
        """EC suggestions + enzyme candidates for a reaction (works with no EC annotated).

        Databases annotate EC numbers patchily, so "which enzyme performs this?" is
        usually unanswerable from the table alone. This looks the EC up from the
        reaction's own cross-references, or infers it from an identical reaction
        elsewhere, then offers the UniProt sequences that could do the job.
        """
        from ..core import enzymes as ez
        db = getattr(self, "_pathway_db", None) or self._combined_included_model()
        host_name = ""
        if self.project is not None and self.project.model is not None:
            host_name = getattr(self.project.model, "name", "") or self.project.model.id

        ok, sugg = run_busy(
            self, f"Looking up EC numbers for {rxn.id}…",
            lambda: ez.suggest_ec_numbers(rxn, db, online=True),
            title="Suggest EC numbers", cancelable=True)
        if not ok:
            if not was_cancelled(sugg):
                QMessageBox.warning(self, "EC lookup failed", str(sugg))
            return
        equation = ""
        try:
            from ..core.network_graph import clean_label
            equation = clean_label(rxn.build_reaction_string(use_metabolite_names=True))
        except Exception:  # noqa: BLE001
            pass
        from .dialogs.enzyme_dialog import EnzymeDialog
        EnzymeDialog(self, rxn, sugg, host_name=host_name, equation=equation).exec()

    def _on_pathway_result_context(self, cells: list, global_pos) -> None:
        """Right-click on a predicted-pathway reaction: show its database details."""
        db = getattr(self, "_pathway_db", None)
        if db is None:
            return
        rid = next((c for c in cells if db.reactions.has_id(c)), None)
        # the table may show a suggested id; fall back to matching by any cell
        if rid is None:
            for c in cells:
                if db.reactions.has_id(c):
                    rid = c
                    break
        if rid is None:
            return
        rxn = db.reactions.get_by_id(rid)
        menu = QMenu(self)
        menu.addAction("Reaction information…", lambda: self._show_reaction_details(rxn))
        # Most database reactions carry no EC, which leaves the user with no route to an
        # enzyme. Offer suggestions (annotated / cross-referenced / inferred) plus the
        # UniProt candidates that could perform the step.
        menu.addAction("Suggest EC numbers && enzymes…",
                       lambda _=False, r=rxn: self._suggest_enzymes(r))
        # Offer H+/H2O balancing when this suggested reaction is unbalanced but fixable.
        from ..core import balancing
        if balancing.plan_water_proton_changes(rxn) is not None:
            menu.addAction("Balance with H⁺/H₂O…",
                           lambda _=False, r=rxn: self._balance_pathway_reaction(r))
        # …and offer to close the gap when the step could not be checked at all.
        if not balancing.assess_reaction(rxn).checkable:
            menu.addAction("Fetch missing formula…",
                           lambda _=False, r=rxn: self._fetch_pathway_formulas(r))
        if self.project is not None and self.project.model.reactions.has_id(rid):
            for steps in (1, 2):
                menu.addAction(f"Show in map ({steps} step{'s' if steps > 1 else ''})",
                               lambda _=False, s=steps, r=rid: self._show_in_map(r, s))
        menu.exec(global_pos)

    def _balance_current_pathway_steps(self) -> None:
        """Balance every fixable unbalanced step in the current pathway result, each
        after user confirmation of the exact H+/H2O change (fix 12)."""
        from ..core import balancing
        db = getattr(self, "_pathway_db", None)
        result = self.pathway_panel.current_result()
        if db is None or result is None or result.reactions.empty:
            return
        fixed, skipped = 0, 0
        for rid in list(result.reactions["reaction"]):
            if not db.reactions.has_id(rid):
                continue
            rxn = db.reactions.get_by_id(rid)
            changes = balancing.plan_water_proton_changes(rxn)
            if not changes:
                continue
            summary = "\n".join(f"  • add {abs(c['delta']):g} {c['metabolite_id']} "
                                f"to the {c['side']} side" for c in changes)
            if QMessageBox.question(
                    self, f"Balance {rid}?",
                    f"Add (aqueous, ~pH 7):\n\n{summary}\n\nApply to {rid}?",
                    QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
                balancing.apply_changes(rxn, changes)
                mask = result.reactions["reaction"] == rid
                result.reactions.loc[mask, "equation"] = pathway_design._equation_with_names(rxn)
                fixed += 1
            else:
                skipped += 1
        self.pathway_panel.refresh_current_result()
        QMessageBox.information(
            self, "Balancing complete",
            f"Balanced {fixed} reaction(s)" + (f", skipped {skipped}" if skipped else "")
            + ". Reactions whose imbalance involves other elements were left unchanged "
            "(they need a missing substrate/cofactor or a corrected formula).")

    def _fetch_pathway_formulas(self, rxn) -> None:
        """Fetch the formulas that stop a suggested step being balance-checked.

        The database gap is the user's to close, but not their job to research: look the
        formula up from BiGG, then from the compound's own cross-references, write it
        into the loaded database, and re-run the balance check so the route's verdict
        changes from "cannot be checked" to a real answer — balanced or not.
        """
        from ..core import balancing, chemistry
        db = getattr(self, "_pathway_db", None)
        model = getattr(rxn, "model", None) or db
        if model is None:
            return
        missing = [m.id for m in balancing.missing_formula_metabolites(rxn)]
        if not missing:
            return
        ok, report = run_busy(
            self, f"Looking up {len(missing)} formula(s)…",
            lambda: chemistry.fetch_missing_formulas(model, missing, online=True),
            title="Fetch missing formula", cancelable=True)
        if not ok:
            if not was_cancelled(report):
                QMessageBox.warning(self, "Could not fetch formulas", str(report))
            return
        result = self.pathway_panel.current_result()
        if result is not None and db is not None and result.reaction_ids:
            self._recheck_balance(result, db, result.reaction_ids)
            self.pathway_panel.refresh_current_result()
        verdict = balancing.assess_reaction(rxn)
        QMessageBox.information(
            self, "Formula lookup",
            report.sentence() + f"\n\n{rxn.id}: {verdict.sentence()}.")

    def _balance_pathway_reaction(self, rxn) -> None:
        """Balance a suggested pathway reaction with H+/H2O after user confirmation
        (fix 12). Edits the reaction in the loaded database so it is balanced when
        the pathway is applied, and refreshes the shown balance status."""
        from ..core import balancing
        changes = balancing.plan_water_proton_changes(rxn)
        if not changes:
            QMessageBox.information(
                self, "Cannot auto-balance",
                "This reaction's imbalance can't be closed with water and protons alone "
                "(other elements, or an inconsistent residual). It likely needs a missing "
                "substrate/cofactor or a corrected formula.")
            return
        summary = "\n".join(f"  • add {abs(c['delta']):g} {c['metabolite_id']} "
                            f"to the {c['side']} side" for c in changes)
        if QMessageBox.question(
                self, "Confirm balancing changes",
                f"Balance {rxn.id} by adding (aqueous, ~pH 7):\n\n{summary}\n\n"
                "This assumes the imbalance is protonation/hydration bookkeeping. "
                "Apply these changes to the suggested reaction?",
                QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return
        balancing.apply_changes(rxn, changes)
        # Refresh the shown equation for this reaction, if it's in the current result.
        result = self.pathway_panel.current_result()
        if result is not None and not result.reactions.empty:
            df = result.reactions
            mask = df["reaction"] == rxn.id
            if mask.any() and "equation" in df.columns:
                df.loc[mask, "equation"] = pathway_design._equation_with_names(rxn)
                self.pathway_panel.refresh_current_result()
        QMessageBox.information(self, "Reaction balanced",
                                f"{rxn.id} is now mass/charge-balanced and will be applied "
                                "balanced when you add the pathway.")

    def _show_reaction_details(self, rxn) -> None:
        """Update the Selection panel and open a summary popup for a reaction."""
        flux = None
        if self._last_result is not None and self._last_result.is_optimal:
            flux = self._last_result.fluxes.get(rxn.id)
        self.info.show_reaction(rxn, flux)
        import urllib.parse
        ec = databases.reaction_ec_numbers(rxn)
        # Each EC number links to a UniProt search for reviewed enzymes in that class.
        ec_links = "—"
        if ec:
            parts = []
            for e in ec:
                q = urllib.parse.quote(f"(ec:{e}) AND (reviewed:true)")
                parts.append(f"<a href='https://www.uniprot.org/uniprotkb?query={q}'>{e}</a>")
            ec_links = ", ".join(parts)
        from ..core.network_graph import clean_label, display_reaction_name
        rows = [
            ("Name", display_reaction_name(rxn) or "—"),
            ("Equation", clean_label(rxn.build_reaction_string(use_metabolite_names=True))),
            ("Bounds", f"[{rxn.lower_bound:g}, {rxn.upper_bound:g}]"),
            ("Direction", "reversible (↔)" if rxn.reversibility else "irreversible (→)"),
            ("Subsystem", rxn.subsystem or "—"),
            ("EC number", ec_links),
            ("Gene rule", rxn.gene_reaction_rule or "—"),
        ]
        if flux is not None:
            rows.append(("Current flux", f"{flux:.6g}"))
        body = "".join(
            f"<tr><td style='padding:2px 14px 2px 0;color:#5f6368'>{k}</td>"
            f"<td style='padding:2px 0'>{v}</td></tr>" for k, v in rows)
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Reaction · {rxn.id}")
        dlg.resize(560, 480)
        lay = QVBoxLayout(dlg)
        label = QLabel(f"<h3 style='margin-top:0'>{rxn.id}</h3><table>{body}</table>")
        label.setTextInteractionFlags(Qt.TextBrowserInteraction)
        label.setOpenExternalLinks(True)
        label.setWordWrap(True)
        label.setMargin(12)
        label.setAlignment(Qt.AlignTop)
        lay.addWidget(label)
        # Graphical depiction of the reaction (structure images = structure images).
        from .widgets.reaction_graphic import ReactionGraphicWidget
        gfx_title = QLabel("Reaction drawing")
        gfx_title.setStyleSheet("color:#5f6368; padding:6px 12px 0;")
        lay.addWidget(gfx_title)
        lay.addWidget(ReactionGraphicWidget(rxn, dlg), 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        dlg.exec()

    def _edit_metabolite(self, met) -> None:
        """Rename a metabolite (id and/or display name), updating every reaction it
        appears in. Useful when a database left it with an opaque id like C20413_c."""
        if self.project is None or met is None:
            return
        from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFormLayout, QLineEdit)
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Edit metabolite · {met.id}")
        form = QFormLayout(dlg)
        comp = met.compartment or (met.id.rsplit("_", 1)[-1] if "_" in met.id else "")
        base = met.id[:-(len(comp) + 1)] if comp and met.id.endswith("_" + comp) else met.id
        id_edit = QLineEdit(base)
        name_edit = QLineEdit(met.name or "")
        form.addRow("Identifier (without compartment):", id_edit)
        form.addRow(f"Compartment:", QLabel(comp or "—"))
        form.addRow("Descriptive name:", name_edit)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, dlg)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        form.addRow(buttons)
        if dlg.exec() != QDialog.Accepted:
            return
        new_base = id_edit.text().strip()
        new_name = name_edit.text().strip()
        import re as _re
        new_base = _re.sub(r"[^A-Za-z0-9_]+", "_", new_base).strip("_")
        if not new_base:
            return
        new_id = f"{new_base}_{comp}" if comp else new_base
        old_id = met.id

        def _apply(m):
            target = m.metabolites.get_by_id(old_id)
            if new_id != old_id and not m.metabolites.has_id(new_id):
                target.id = new_id       # cobra updates every reaction referencing it
            target.name = new_name or target.name
            m.repair()
        if new_id != old_id and self.project.model.metabolites.has_id(new_id):
            QMessageBox.warning(self, "Edit metabolite",
                                f"A metabolite '{new_id}' already exists.")
            return
        self.project.apply_edit(_apply)
        self._refresh_after_edit(f"Renamed metabolite {old_id} → {new_id}.")

    def _show_metabolite_details(self, met) -> None:
        """Update the Selection panel and open a summary popup for a metabolite."""
        self.info.show_metabolite(met)
        ec_rows = ""
        for key in ("kegg.compound", "bigg.metabolite", "chebi", "metanetx.chemical"):
            val = met.annotation.get(key) if isinstance(met.annotation, dict) else None
            if val:
                disp = ", ".join(val) if isinstance(val, (list, tuple)) else str(val)
                ec_rows += (f"<tr><td style='padding:2px 14px 2px 0;color:#5f6368'>{key}</td>"
                            f"<td style='padding:2px 0'>{disp}</td></tr>")
        rxns = ", ".join(r.id for r in met.reactions) or "—"
        rows = [
            ("Name", met.name or "—"),
            ("Formula", met.formula or "—"),
            ("Charge", str(met.charge) if met.charge is not None else "—"),
            ("Compartment", met.compartment or "—"),
            ("# Reactions", str(len(met.reactions))),
        ]
        body = "".join(
            f"<tr><td style='padding:2px 14px 2px 0;color:#5f6368'>{k}</td>"
            f"<td style='padding:2px 0'>{v}</td></tr>" for k, v in rows)
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Metabolite · {met.id}")
        dlg.resize(460, 360)
        lay = QVBoxLayout(dlg)
        scroll = QScrollArea(dlg)
        scroll.setWidgetResizable(True)
        label = QLabel(
            f"<h3 style='margin-top:0'>{met.id}</h3><table>{body}{ec_rows}</table>"
            f"<p style='color:#5f6368'>Reactions:</p><p>{rxns}</p>")
        label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        label.setWordWrap(True)
        label.setMargin(12)
        label.setAlignment(Qt.AlignTop)
        scroll.setWidget(label)
        lay.addWidget(scroll)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        dlg.exec()

    def _show_model_summary_on_load(self) -> None:
        """A compact Model Summary popup shown right after a model loads, since the
        Information panel is hidden by default. Reminds the user where to find it."""
        if self.project is None or self.project.model is None:
            return
        try:
            summary = io_models.summarize(self.project.model)
            rows = "".join(
                f"<tr><td style='padding:2px 16px 2px 0;color:#5f6368'>{label}</td>"
                f"<td style='padding:2px 0'>{value}</td></tr>"
                for label, value in summary.as_rows())
        except Exception:  # noqa: BLE001 - a summary must never block loading
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("Model summary")
        dlg.resize(440, 360)
        lay = QVBoxLayout(dlg)
        scroll = QScrollArea(dlg)
        scroll.setWidgetResizable(True)
        body = QLabel(
            f"<h3 style='margin-top:0'>Model loaded</h3><table>{rows}</table>"
            "<p style='color:#5f6368;margin-top:10px'>The <b>Information</b> panel is "
            "hidden to keep the window tidy — reopen it any time from "
            "<b>View ▸ Information</b>. Selecting a reaction or metabolite fills it.</p>")
        body.setTextInteractionFlags(Qt.TextSelectableByMouse)
        body.setWordWrap(True)
        body.setMargin(12)
        body.setAlignment(Qt.AlignTop)
        scroll.setWidget(body)
        lay.addWidget(scroll)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        dlg.exec()

    def _show_model_info(self) -> None:
        if self.project is None or self.project.model is None:
            QMessageBox.information(self, "Model info", "Open a model first.")
            return
        summary = io_models.summarize(self.project.model)
        rows = "".join(
            f"<tr><td style='padding:2px 14px 2px 0;color:#5f6368'>{label}</td>"
            f"<td style='padding:2px 0'>{value}</td></tr>"
            for label, value in summary.as_rows()
        )
        dlg = QDialog(self)
        dlg.setWindowTitle("Model info")
        dlg.resize(460, 360)
        lay = QVBoxLayout(dlg)
        scroll = QScrollArea(dlg)
        scroll.setWidgetResizable(True)
        body = QLabel(f"<h3 style='margin-top:0'>General information</h3><table>{rows}</table>")
        body.setTextInteractionFlags(Qt.TextSelectableByMouse)
        body.setWordWrap(True)
        body.setMargin(12)
        body.setAlignment(Qt.AlignTop)
        scroll.setWidget(body)
        lay.addWidget(scroll)
        buttons = QDialogButtonBox(QDialogButtonBox.Close, dlg)
        buttons.rejected.connect(dlg.reject)
        buttons.accepted.connect(dlg.accept)
        lay.addWidget(buttons)
        dlg.exec()

    def _about(self) -> None:
        QMessageBox.about(
            self, "About GSM ToolBox",
            f"<b>{__app_name__}</b> v{__version__}<br><br>"
            "A user-friendly desktop app for genome-scale metabolic modeling "
            "and metabolic engineering.<br><br>Built on COBRApy.")

    def _open_preferences(self) -> None:
        """Settings ▸ Preferences. Optional suites (MDF) are toggled here; the pathway
        panel is refreshed afterwards so newly-enabled controls appear immediately."""
        from .dialogs.preferences_dialog import PreferencesDialog
        from ..core import preferences
        if PreferencesDialog(self).exec():
            preferences.reload()
            # MDF controls are gated on the preference: re-evaluate button visibility.
            self.pathway_panel.set_mdf_enabled(self._mdf_enabled())
            self.pathway_panel.refresh_action_visibility()
            # The quick-access bar and the look of the window are both preferences, so
            # they are re-applied here rather than waiting for a restart.
            self._build_toolbar()
            self._update_action_availability()
            self._apply_appearance()

    @staticmethod
    def _mdf_enabled() -> bool:
        from ..core import thermodynamics as thermo
        return thermo.mdf_suite_enabled()

    def _check_for_updates(self) -> None:
        """Read-only report of the toolbox, its dependencies, and its data assets, with
        suggested manual update commands. Never installs anything itself."""
        from .dialogs.update_dialog import UpdateDialog
        UpdateDialog(self).exec()
