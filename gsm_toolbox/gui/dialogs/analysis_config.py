"""Declarative per-analysis configuration dialogs.

Each analysis declares a list of :class:`Param` descriptors; :class:`AnalysisConfigDialog`
renders the matching widgets (reaction/metabolite pickers, spin boxes, check boxes,
choices) and returns the chosen values. This gives every analysis a settings pop-up
without hand-coding a dialog per method.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from ...core import editing


@dataclass
class Param:
    key: str
    label: str
    kind: str            # reaction | metabolite | int | float | bool | choice
    default: Any = None
    minimum: float = -1e9
    maximum: float = 1e9
    decimals: int = 3
    step: float = 1.0
    choices: List[str] = field(default_factory=list)
    help: str = ""


class AnalysisConfigDialog(QDialog):
    def __init__(self, title: str, params: List[Param], model: cobra.Model, parent=None,
                 description: str = ""):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(460)
        self._params = params
        self._widgets: Dict[str, QWidget] = {}

        rxn_ids = [r.id for r in model.reactions]
        met_ids = [m.id for m in model.metabolites]

        desc_label = None
        if description:
            desc_label = QLabel(description)
            desc_label.setWordWrap(True)
            desc_label.setTextFormat(Qt.RichText)
            desc_label.setStyleSheet(
                "background:#EEF3FB; border:1px solid #D6E0F0; border-radius:6px; "
                "padding:8px; color:#33415c;")

        form_host = QWidget()
        form = QFormLayout(form_host)
        for p in params:
            w = self._make_widget(p, rxn_ids, met_ids)
            self._widgets[p.key] = w
            label = QLabel(p.label)
            if p.help:
                label.setToolTip(p.help)
                w.setToolTip(p.help)
            form.addRow(label, w)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        scroll.setWidget(form_host)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        if desc_label is not None:
            layout.addWidget(desc_label)
        intro = QLabel("Adjust the settings for this analysis, then click OK.")
        intro.setStyleSheet("color: #5f6368;")
        layout.addWidget(scroll, 1)
        layout.addWidget(intro)
        layout.addWidget(buttons)
        if len(params) > 7:
            self.setMinimumHeight(560)

    def _make_widget(self, p: Param, rxn_ids, met_ids) -> QWidget:
        if p.kind in ("reaction", "metabolite", "choice"):
            combo = QComboBox()
            items = p.choices if p.kind == "choice" else (rxn_ids if p.kind == "reaction" else met_ids)
            combo.addItems([str(i) for i in items])
            if p.kind == "choice":
                combo.setEditable(False)
            else:
                # Editable search box: contains-filter popup + select-all on focus,
                # so typing replaces/filters instead of appending (#2).
                from ..widgets.dialog_util import configure_search_combo
                configure_search_combo(combo)
            if p.default is not None and str(p.default) in [str(i) for i in items]:
                combo.setCurrentText(str(p.default))
            elif p.default is not None and p.kind != "choice":
                combo.setCurrentText(str(p.default))
            return combo
        if p.kind == "int":
            spin = QSpinBox()
            spin.setRange(int(p.minimum), int(p.maximum))
            spin.setValue(int(p.default if p.default is not None else 0))
            return spin
        if p.kind == "float":
            spin = QDoubleSpinBox()
            spin.setDecimals(p.decimals)
            spin.setRange(float(p.minimum), float(p.maximum))
            spin.setSingleStep(p.step)
            spin.setValue(float(p.default if p.default is not None else 0.0))
            return spin
        if p.kind == "bool":
            chk = QCheckBox()
            chk.setChecked(bool(p.default))
            return chk
        raise ValueError(f"Unknown param kind: {p.kind}")

    def values(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        for p in self._params:
            w = self._widgets[p.key]
            if isinstance(w, QComboBox):
                out[p.key] = w.currentText().strip()
            elif isinstance(w, QSpinBox):
                out[p.key] = w.value()
            elif isinstance(w, QDoubleSpinBox):
                out[p.key] = w.value()
            elif isinstance(w, QCheckBox):
                out[p.key] = w.isChecked()
        return out


# ----------------------------------------------------------------------
# Per-analysis parameter schemas
# ----------------------------------------------------------------------
def _first_exchange(model: cobra.Model) -> Optional[str]:
    for r in model.reactions:
        if r.id.startswith("EX_"):
            return r.id
    return model.reactions[0].id if model.reactions else None


# Plain-language explanations shown at the top of an analysis's settings dialog.
ANALYSIS_DESCRIPTIONS: Dict[str, str] = {
    "overproduction": (
        "<b>Metabolite overproduction.</b> You choose a <i>metabolite</i> to accumulate. "
        "The tool temporarily adds a demand reaction that drains it, which forces the network "
        "to <i>net-produce</i> that metabolite — this works even when the compound is a "
        "reactant or is made by reversible reactions. It then runs an <b>FSEOF</b> scan: the "
        "enforced demand is increased step by step while growth is maximised, and every "
        "reaction's flux is tracked.<br><br>"
        "<b>Reading the results:</b> reactions whose flux rises as production is pushed are "
        "<i>overexpression</i> targets (amplify them); reactions whose flux falls toward zero "
        "are <i>knockdown/deletion</i> targets. The <b>direction</b> column tells you which "
        "way each reaction must run (forward/reverse). Targets are ranked by how much their "
        "flux changes."),
    "knockout": (
        "<b>Knockout strain design.</b> Searches for reaction knockouts that couple production "
        "of the target to growth. OptKnock/RobustKnock/OptCouple are exact methods; the "
        "heuristic is a faster evolutionary search. Each result lists the knockout set with "
        "the predicted growth and product flux."),
    "fseof": (
        "<b>FSEOF.</b> Enforces an increasing fraction of the maximum flux through the target "
        "<i>reaction</i> while maximising growth, then classifies each reaction as an "
        "overexpression or knock-down target from how its flux responds."),
}


def description_for(analysis_id: str) -> str:
    return ANALYSIS_DESCRIPTIONS.get(analysis_id, "")


# Analyses that already define their own transport/exchange controls, or where a
# reaction-exclusion toggle is meaningless, are skipped when auto-appending them.
_SKIP_EXCLUSION_TOGGLE = {"knockout"}


def params_for(analysis_id: str, model: cobra.Model, context: dict) -> List[Param]:
    """Return the configurable parameters for an analysis.

    The two transport/exchange exclusion toggles are appended to every analysis
    (defaulting to excluded) so a settings dialog always appears with them.
    """
    params = list(_params_for_base(analysis_id, model, context))
    if analysis_id not in _SKIP_EXCLUSION_TOGGLE:
        params.append(Param(
            "exclude_transport", "Exclude transport reactions from results", "bool", True,
            help="Hide membrane-transport reactions from this analysis's results / targets."))
        params.append(Param(
            "exclude_exchange", "Exclude exchange reactions from results", "bool", True,
            help="Hide exchange / boundary reactions from this analysis's results / targets."))
    return params


def _params_for_base(analysis_id: str, model: cobra.Model, context: dict) -> List[Param]:
    """Return the analysis-specific parameters (empty = none)."""
    biomass = editing.guess_biomass_reaction(model) or ""
    product = _first_exchange(model)
    selected = context.get("selected_reaction") or product

    if analysis_id == "fva":
        return [
            Param("fraction_of_optimum", "Fraction of optimum to maintain (0–1)", "float",
                  1.0, 0.0, 1.0, 2, 0.05,
                  help="Constrain the objective to at least this fraction of its optimum."),
            Param("loopless", "Loopless (remove thermodynamically infeasible cycles)", "bool", False),
        ]
    if analysis_id == "production_envelope":
        return [
            Param("target", "Target product reaction", "reaction", product),
            Param("points", "Number of points", "int", 20, 5, 200),
        ]
    if analysis_id == "robustness":
        return [
            Param("control", "Reaction to vary", "reaction", selected),
            Param("auto_range", "Auto range from flux variability (recommended)", "bool", True,
                  help="Scan across the reaction's feasible flux range."),
            Param("lower", "Lower flux (if not auto)", "float", 0.0, -1e6, 1e6, 3),
            Param("upper", "Upper flux (if not auto)", "float", 10.0, -1e6, 1e6, 3),
            Param("points", "Number of points", "int", 25, 5, 200),
        ]
    if analysis_id == "phase_plane":
        return [
            Param("reaction_x", "First reaction (x axis)", "reaction", selected),
            Param("reaction_y", "Second reaction (y axis)", "reaction", product),
            Param("points", "Points per axis", "int", 15, 4, 40),
        ]
    if analysis_id == "quality_report":
        return [Param("include_blocked", "Find blocked reactions (slower)", "bool", True)]
    if analysis_id == "knockout":
        from ...core.analysis import strain_design as sd
        solver_choices = ["Auto"] + sd.available_solvers()
        return [
            Param("method", "Method", "choice", "OptKnock",
                  choices=["OptKnock", "RobustKnock", "OptCouple", "Heuristic (evolutionary)"],
                  help="OptKnock/RobustKnock/OptCouple are exact MILP methods. Heuristic is a "
                       "fast evolutionary search for when exact methods are intractable."),
            Param("product", "Target product reaction (to maximize)", "reaction", product,
                  help="The reaction whose flux is coupled to growth — choose the EXCHANGE "
                       "reaction that secretes your product (e.g. EX_btoh_e), not an internal "
                       "step. Maximizing an internal reaction can be misleading."),
            Param("biomass", "Biomass / growth reaction", "reaction", biomass),
            Param("max_knockouts", "Maximum number of knockouts", "int", 3, 1, 15,
                  help="Keep small (≤3–4); the search space grows exponentially."),
            Param("max_solutions", "Maximum solutions to return", "int", 5, 1, 50),
            Param("min_growth_fraction", "Minimum growth (fraction of wild-type)", "float",
                  0.1, 0.0, 1.0, 2, 0.05,
                  help="Forces viable designs. If no solutions are found, lower it (e.g. 0.01)."),
            Param("solver", "Solver (exact methods)", "choice", "Auto", choices=solver_choices,
                  help="SCIP is fast for OptKnock; RobustKnock uses GLPK. Ignored by the heuristic."),
            Param("approach", "Search approach (exact methods)", "choice", "Best (optimal)",
                  choices=["Best (optimal)", "Any (faster)", "Diverse set"],
                  help="'Any' returns a quick suboptimal design; 'Diverse set' returns varied designs."),
            Param("exclude_exchanges", "Protect exchange reactions from knockout", "bool", True),
            Param("exclude_transport", "Protect transport reactions from knockout", "bool", False),
            Param("exclude_blocked", "Ignore blocked / no-flux reactions", "bool", True),
            Param("time_limit", "Time limit (seconds)", "int", 120, 10, 3600),
            Param("population", "Population size (heuristic only)", "int", 40, 10, 300),
            Param("generations", "Generations (heuristic only)", "int", 25, 5, 300),
        ]
    if analysis_id == "flux_sampling":
        return [
            Param("n_samples", "Number of samples", "int", 500, 50, 5000,
                  help="More samples give smoother distributions but take longer. "
                       "Restrict to a category (Scope) to focus the reactions shown."),
        ]
    if analysis_id == "fseof":
        return [
            Param("product", "Target product reaction", "reaction", product,
                  help="The reaction whose flux is progressively enforced — choose the EXCHANGE "
                       "reaction that secretes your product (e.g. EX_btoh_e). FSEOF then finds "
                       "reactions to amplify/knock down to push flux toward it."),
            Param("n_steps", "Number of enforced-flux steps", "int", 10, 3, 50),
            Param("tolerance", "Monotonicity tolerance (per-step slack)", "float",
                  1e-6, 0.0, 1.0, 8, 1e-6,
                  help="Slack allowed when testing strict monotonicity — raise it to keep "
                       "near-flat real trends."),
            Param("include_trending", "Also report strong non-monotone trends", "bool", True,
                  help="Include reactions correlated with the enforced flux even if not "
                       "strictly monotone (marked '(trend)')."),
        ]
    if analysis_id == "overproduction":
        return [
            Param("metabolite", "Metabolite to accumulate / overproduce", "metabolite", None,
                  help="Search for the target metabolite. A demand reaction is added to "
                       "force its net production, then reactions to amplify/knock down are found."),
            Param("n_steps", "Number of enforced-flux steps", "int", 10, 3, 50),
            Param("tolerance", "Monotonicity tolerance (per-step slack)", "float",
                  1e-6, 0.0, 1.0, 8, 1e-6),
            Param("include_trending", "Also report strong non-monotone trends", "bool", True),
            Param("secretion", "Enforce via real secretion (transport + exchange)", "bool", False,
                  help="Instead of a bare DM_ demand drain, add a transport step to the "
                       "extracellular compartment and an exchange reaction, so the target "
                       "reflects realistic export."),
            Param("gas", "Volatile product (add a gas-phase sink)", "bool", False,
                  help="For products lost to the gas phase."),
        ]
    if analysis_id == "gimme":
        median = context.get("expression_median", 1.0)
        return [
            Param("threshold", "Expression threshold (below = penalized)", "float",
                  round(float(median), 3), 0.0, 1e9, 3, 1.0),
            Param("objective_fraction", "Fraction of optimum to maintain (0–1)", "float",
                  0.9, 0.0, 1.0, 2, 0.05),
        ]
    if analysis_id == "atpm_sensitivity":
        atpm = "ATPM" if model.reactions.has_id("ATPM") else ""
        return [
            Param("atpm_id", "ATP-maintenance reaction", "reaction", atpm),
            Param("points", "Number of points", "int", 20, 5, 100),
        ]
    if analysis_id == "mutant":
        return [Param("method", "Method", "choice", "MOMA (minimal adjustment)",
                      choices=["MOMA (minimal adjustment)", "ROOM (regulatory on/off)"],
                      help="MOMA: flux state closest to wild-type. ROOM: fewest reactions "
                           "changed. Knockouts are the reactions selected in the Explorer.")]
    if analysis_id == "efm":
        return [Param("max_modes", "Maximum modes to list", "int", 200, 10, 5000)]
    if analysis_id in ("gapfill_growth", "gapfill_metabolite"):
        params = []
        if analysis_id == "gapfill_metabolite":
            params.append(Param("metabolite", "Metabolite to produce", "metabolite", None))
        params.append(Param("lower_bound", "Minimum flux to restore", "float", 0.05, 0.0, 1000.0, 3))
        return params
    return []
