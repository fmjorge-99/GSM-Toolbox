"""The central 'Analysis' tab: pick an analysis, optionally restrict to a category,
and view results as a table or plot.

The panel only declares *what* the user wants to run (via ``run_requested``); the
main window gathers any needed parameters and runs the job on a worker thread.
"""

from __future__ import annotations

from typing import List, Optional

import pandas as pd
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..views.plot_view import PlotView
from ..views.results_view import ResultsView
from ..widgets.wrap_button import WrapButton
from .objective_bar import ObjectiveBar


# (analysis_id, button label, tooltip explaining the method in plain language)
SIMULATION = [
    ("fba", "Flux Balance Analysis (FBA)",
     "Predict the flux through every reaction that maximizes the objective (e.g. growth)."),
    ("pfba", "Parsimonious FBA (pFBA)",
     "Like FBA, but among optimal solutions picks the one using the least total enzyme flux."),
    ("fva", "Flux Variability Analysis (FVA)",
     "For each reaction, the min and max flux possible while staying near the optimum."),
    ("shadow_prices", "Shadow prices",
     "Which metabolites limit the objective — how much the objective would change per unit of each."),
    ("mutant", "Mutant phenotype (MOMA/ROOM)…",
     "Predict a knockout mutant's fluxes (often more realistic than FBA). Select the "
     "reactions to knock out in the Explorer first."),
]
ESSENTIALITY = [
    ("single_reaction_deletion", "Single reaction deletion",
     "Knock out each reaction one at a time and report the effect on growth (finds essential reactions)."),
    ("single_gene_deletion", "Single gene deletion",
     "Knock out each gene one at a time and report the effect on growth (finds essential genes)."),
]
PHENOTYPE = [
    ("production_envelope", "Production envelope",
     "Feasible range of a product's flux vs growth — shows growth/production trade-offs."),
    ("robustness", "Robustness analysis",
     "Fix one reaction across a range and watch the objective — reveals tipping points."),
    ("phase_plane", "Phenotypic phase plane",
     "Vary two reactions and map the objective surface (e.g. carbon vs oxygen uptake)."),
    ("flux_sampling", "Flux sampling (distributions)",
     "Sample the feasible flux space and show per-reaction violin distributions — "
     "conveys uncertainty and alternate optima."),
]
STRAIN_DESIGN = [
    ("overproduction", "Overproduce a metabolite…",
     "Pick a metabolite to accumulate; the toolbox finds which reactions to amplify or "
     "knock down (and in which direction) to maximize it — handles reversible reactions and "
     "targets where the metabolite is a reactant."),
    ("knockout", "Knockout strain design…",
     "Design reaction knockouts that couple a product to growth. Choose the method "
     "(OptKnock / RobustKnock / OptCouple / Heuristic), solver and parameters in the dialog."),
    ("fseof", "FSEOF (amplification targets)",
     "Scan increasing enforced product flux to find over-expression and knock-down targets."),
]
OMICS = [
    ("prepare_omics", "Prepare omics dataset…",
     "Turn a raw transcriptomics / proteomics / metabolomics table (any format) into the "
     "(id, value) form the model needs, mapping its identifiers onto the model's namespace."),
    ("load_expression", "Load expression data…",
     "Load a ready two-column gene-expression table (CSV/TSV: gene id, value) for eFlux/GIMME."),
    ("eflux", "eFlux (expression-constrained)",
     "Scale reaction bounds by gene expression, then simulate — a context-specific flux state."),
    ("gimme", "GIMME (context-specific)…",
     "Minimize flux through low-expression reactions while meeting the objective."),
    ("atpm_sensitivity", "ATP maintenance sensitivity",
     "Scan the non-growth ATP maintenance demand and plot its effect on growth."),
]
PATHWAYS = [
    ("efm", "Elementary flux modes (category)",
     "Enumerate the minimal steady-state flux routes through a small category/subnetwork."),
    ("community_growth", "Community member growth",
     "Run FBA on a community model and report each organism's growth (build a community first)."),
]
CURATION = [
    ("quality_report", "Quality report",
     "Check the model for mass/charge imbalances, missing gene rules, dead-end metabolites and blocked reactions."),
    ("blocked_reactions", "Find blocked reactions",
     "List reactions that cannot carry any flux under the current constraints."),
    ("gapfill_growth", "Gap-fill for growth…",
     "Suggest reactions from a universal model that would let a non-growing model grow."),
    ("gapfill_metabolite", "Gap-fill for a metabolite…",
     "Suggest reactions that would let the model produce a chosen metabolite."),
]


class AnalysisPanel(QWidget):
    run_requested = Signal(str)  # analysis id
    save_requested = Signal()    # save the currently-shown result
    display_fluxes_requested = Signal()   # draw the current flux result on the network
    visualize_requested = Signal()        # open the Plot Gallery for the current result

    def __init__(self):
        super().__init__()

        # Objective settings live in a popup opened from the top (#C13), not a
        # left-column panel — keeping the analysis controls compact.
        self.objective_bar = ObjectiveBar()

        # Category scope selector — at the TOP of the tools (#T4).
        scope_box = QGroupBox("Scope")
        scope_layout = QHBoxLayout(scope_box)
        scope_layout.addWidget(QLabel("Run on:"))
        self.scope_combo = QComboBox()
        self.scope_combo.addItem("Whole model", userData=None)
        self.scope_combo.setToolTip(
            "Restrict the analysis to a category. Connections to the rest of the network "
            "are treated as free input/output fluxes."
        )
        scope_layout.addWidget(self.scope_combo, 1)

        controls = QVBoxLayout()
        controls.addWidget(scope_box)
        controls.addWidget(self._make_group("Simulation", SIMULATION, primary_first=True))
        controls.addWidget(self._make_group("Strain design", STRAIN_DESIGN))
        # Phenotype & Essentiality combined into one category (#C14).
        controls.addWidget(self._make_group("Phenotype & essentiality",
                                            PHENOTYPE + ESSENTIALITY))
        controls.addWidget(self._make_group("Omics & energy", OMICS))
        controls.addWidget(self._make_group("Pathways & community", PATHWAYS))
        controls.addWidget(self._make_group("Model quality & curation", CURATION))
        controls.addStretch(1)

        controls_inner = QWidget()
        controls_inner.setLayout(controls)
        controls_widget = QScrollArea()
        controls_widget.setWidgetResizable(True)
        controls_widget.setWidget(controls_inner)
        controls_widget.setMinimumWidth(220)
        controls_widget.setFrameShape(QScrollArea.NoFrame)

        # Active-physiology banner (Issue 6): objective + carbon/energy source, so a
        # misconfigured medium (e.g. dark heterotrophy on a phototroph) is visible.
        self.physiology_banner = QLabel()
        self.physiology_banner.setWordWrap(True)
        self.physiology_banner.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.physiology_banner.setVisible(False)

        # Results area: Table + Plot tabs
        self.results_view = ResultsView()
        self.plot_view = PlotView()
        self.result_tabs = QTabWidget()
        self.result_tabs.addTab(self.results_view, "Table")
        self.result_tabs.addTab(self.plot_view, "Plot")

        self.flux_map_btn = QPushButton("Display fluxes in network")
        self.flux_map_btn.setToolTip("Draw the reaction fluxes from this result on a network map "
                                     "(width/colour by flux; labels show values and any bounds).")
        self.flux_map_btn.clicked.connect(self.display_fluxes_requested)
        self.flux_map_btn.setVisible(False)
        self.visualize_btn = QPushButton("Visualize…")
        self.visualize_btn.setToolTip("Open the Plot Gallery: publication-grade figures for this "
                                      "result (FVA tornado, FSEOF scan, exchange bars, design "
                                      "comparison, strategy heatmap/waterfall) with SVG/PDF/CSV export.")
        self.visualize_btn.clicked.connect(self.visualize_requested)
        self.save_btn = QPushButton("Save Results…")
        self.save_btn.setToolTip("Save the results shown here to your results folder — a table is "
                                 "saved as CSV, a plot as an image (PNG) and PDF.")
        self.save_btn.clicked.connect(self.save_requested)
        save_row = QHBoxLayout()
        save_row.addStretch(1)
        save_row.addWidget(self.flux_map_btn)
        save_row.addWidget(self.visualize_btn)
        save_row.addWidget(self.save_btn)

        # Top row: active-setup banner + a compact Objective Settings button (#C13).
        self.objective_btn = QPushButton("Objective Settings…")
        self.objective_btn.setToolTip("Choose what to optimise (maximize growth / a product / "
                                      "a custom weighted objective).")
        self.objective_btn.clicked.connect(self._open_objective_settings)
        top_row = QHBoxLayout()
        top_row.addWidget(self.physiology_banner, 1)
        top_row.addWidget(self.objective_btn, 0, Qt.AlignTop)

        results_area = QWidget()
        results_layout = QVBoxLayout(results_area)
        results_layout.setContentsMargins(0, 0, 0, 0)
        results_layout.addLayout(top_row)
        results_layout.addWidget(self.result_tabs, 1)
        results_layout.addLayout(save_row)

        # Draggable divider so the user can rebalance controls vs. results
        # (the results panel used to grow so wide it hid the buttons).
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(controls_widget)
        splitter.addWidget(results_area)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setCollapsible(0, False)
        splitter.setChildrenCollapsible(False)
        splitter.setSizes([300, 700])

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)

    def _make_group(self, title: str, items, primary_first: bool = False) -> QGroupBox:
        box = QGroupBox(title)
        grid = QGridLayout(box)
        for i, (analysis_id, label, tip) in enumerate(items):
            btn = WrapButton(label)     # wraps to width; never crops (#T3)
            btn.setToolTip(tip)
            if primary_first and i == 0:
                btn.setObjectName("primary")
            btn.clicked.connect(lambda _=False, a=analysis_id: self.run_requested.emit(a))
            grid.addWidget(btn, i, 0)
        return box

    # ----- scope -------------------------------------------------------
    def set_categories(self, names: List[str]) -> None:
        current = self.current_category()
        self.scope_combo.blockSignals(True)
        self.scope_combo.clear()
        self.scope_combo.addItem("Whole model", userData=None)
        for name in names:
            self.scope_combo.addItem(f"Category: {name}", userData=name)
        idx = self.scope_combo.findData(current)
        self.scope_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.scope_combo.blockSignals(False)

    def current_category(self) -> Optional[str]:
        return self.scope_combo.currentData()

    def _open_objective_settings(self) -> None:
        """Show the objective editor in a popup (#C13). The ObjectiveBar widget is
        kept live and simply re-parented into a reusable dialog."""
        if getattr(self, "_objective_dialog", None) is None:
            dlg = QDialog(self)
            dlg.setWindowTitle("Objective settings")
            dlg.setMinimumWidth(460)
            lay = QVBoxLayout(dlg)
            lay.addWidget(self.objective_bar)
            bb = QDialogButtonBox(QDialogButtonBox.Close)
            bb.rejected.connect(dlg.reject)
            bb.accepted.connect(dlg.accept)
            lay.addWidget(bb)
            self._objective_dialog = dlg
        self._objective_dialog.show()
        self._objective_dialog.raise_()

    def set_physiology(self, summary) -> None:
        """Update the active-physiology banner (Issue 6). ``summary`` is a
        core.physiology.PhysiologySummary (or None to hide the banner)."""
        if summary is None:
            self.physiology_banner.setVisible(False)
            return
        base = (f"<b>Active setup</b> — objective: <b>{summary.objective}</b> · "
                f"energy: <b>{summary.energy_source}</b> · "
                f"carbon: <b>{summary.carbon_source}</b>")
        if summary.warnings:
            body = base + "<br>⚠ " + "<br>⚠ ".join(summary.warnings)
            style = ("background:#FEF7E0; border:1px solid #F0C36D; border-radius:6px; "
                     "padding:6px 8px; color:#7A5B00;")
        else:
            body = base
            style = ("background:#E8F0FE; border:1px solid #C6D6F5; border-radius:6px; "
                     "padding:6px 8px; color:#1F3D7A;")
        self.physiology_banner.setText(body)
        self.physiology_banner.setStyleSheet(style)
        self.physiology_banner.setVisible(True)

    def set_flux_available(self, available: bool) -> None:
        """Show the 'Display fluxes in network' button when the current result has
        per-reaction flux data."""
        self.flux_map_btn.setVisible(available)

    # ----- results -----------------------------------------------------
    def show_table(self, df: pd.DataFrame, header: str) -> None:
        self.results_view.show_dataframe(df, header)
        self.result_tabs.setCurrentWidget(self.results_view)

    def active_result(self):
        """Return ``(kind, payload, title)`` for the currently-shown result.

        ``kind`` is ``"table"`` (payload = DataFrame), ``"plot"`` (payload =
        PlotView) or ``None`` when there is nothing to save.
        """
        if self.result_tabs.currentWidget() is self.plot_view and self.plot_view.has_content():
            return "plot", self.plot_view, self.plot_view.title
        df = self.results_view.current_dataframe()
        if df is not None and not getattr(df, "empty", True):
            return "table", df, self.results_view.current_header()
        # Fall back to whichever tab has content.
        if self.plot_view.has_content():
            return "plot", self.plot_view, self.plot_view.title
        return None, None, ""
