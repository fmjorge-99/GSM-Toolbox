"""The full feasibility report for a designed pathway (VI.10).

The results table makes it easy to rank routes by predicted flux — the one number that
most often points at the wrong route. This dialog gathers everything the toolbox knows
about a route in one place: the chemistry (skeletal-isomer check, balance), the
thermodynamics, the competing reactions, the flux, and what blocks it if anything does.
"""
from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout, QWidget)

from ...core import feasibility as fz
from .. import style


def _section(title: str) -> QLabel:
    lbl = QLabel(f"<b>{title}</b>")
    lbl.setTextFormat(Qt.RichText)
    lbl.setStyleSheet("margin-top:10px;")
    return lbl


def _body(html: str, colour: str = "") -> QLabel:
    lbl = QLabel(html)
    lbl.setWordWrap(True)
    lbl.setTextFormat(Qt.RichText)
    if colour:
        lbl.setStyleSheet(f"color:{colour};")
    return lbl


class FeasibilityDialog(QDialog):
    """Read-only report; the verdict sentence is repeated at the top."""

    fetch_missing_requested = None    # set by the caller when a gap-fill is offered

    def __init__(self, parent, report: fz.FeasibilityReport, *,
                 on_fetch_missing=None, blockage=None):
        super().__init__(parent)
        self.setWindowTitle(f"Feasibility — {report.target}")
        self.resize(720, 640)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
        self._on_fetch = on_fetch_missing

        outer = QVBoxLayout(self)
        head = QLabel(f"<span style='color:{report.colour}; font-size:13px'>"
                      f"<b>{report.label}</b></span><br>{report.sentence()}")
        head.setWordWrap(True)
        head.setTextFormat(Qt.RichText)
        head.setStyleSheet("padding:8px; background:#F1F3F4; border-radius:4px;")
        outer.addWidget(head)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        v = QVBoxLayout(inner)

        # ---- chemistry ---------------------------------------------------------
        v.addWidget(_section("Chemistry"))
        # Missing structures and formulas are fetched automatically before the checks run,
        # so this reports what was recovered rather than asking the user to go and get it.
        if getattr(report, "backfill_note", ""):
            v.addWidget(_body(f"🔎 {report.backfill_note}", style.TEXT_MUTED))
        if report.isomer_warnings:
            v.addWidget(_body(
                "🚩 <b>A step changes the carbon backbone in a way no single reaction "
                "can.</b> This usually means a database entry is mislabelled — the route "
                "may in fact make a different isomer of your target.", "#C5221F"))
            for w in report.isomer_warnings:
                v.addWidget(_body(f"&nbsp;&nbsp;• {w}"))
        elif report.isomer_checked:
            v.addWidget(_body(
                f"✓ Skeletal-isomer check passed ({report.isomer_coverage} step(s) had "
                "structures to check). No step rearranges the carbon backbone illegally.",
                "#188038"))
        else:
            v.addWidget(_body(
                "The skeletal-isomer check could not run: no step had structures on both "
                "sides. Structures can be fetched from the reaction's cross-references — "
                "until then this route is <b>unchecked</b>, not verified.", style.TEXT_MUTED))

        if not report.balanced:
            v.addWidget(_body("⚠ Some step(s) are mass/charge-<b>unbalanced</b>. Use "
                              "“Balance H⁺/H₂O” or check them by hand.", "#E8710A"))
        else:
            extra = (f" ({len(report.unverified_steps)} step(s) could not be checked "
                     "because a participant has no formula — that is a database gap, not "
                     "an imbalance)" if report.unverified_steps else "")
            v.addWidget(_body(f"✓ Mass and charge balanced{extra}.", "#188038"))
        # Name the participants that defeated the check: "unverified" is only actionable
        # once the user knows which compound to go and get a formula for.
        for rid, labels in list(report.unverified_reasons.items())[:6]:
            if labels:
                v.addWidget(_body(f"&nbsp;&nbsp;• <b>{rid}</b> cannot be checked — "
                                  + ", ".join(labels[:3]) + " has/have no formula.",
                                  style.TEXT_MUTED))

        # ---- thermodynamics ----------------------------------------------------
        v.addWidget(_section("Thermodynamics"))
        if report.mdf is None:
            v.addWidget(_body(
                "No thermodynamic data is available for this route. Either the "
                "thermodynamics suite is disabled (Settings ▸ Preferences), or none of "
                "its reactions could be assigned a ΔrG′°.", style.TEXT_MUTED))
        elif report.mdf_single_reaction:
            v.addWidget(_body(
                "This route is a <b>single reaction</b>, so a Max-min Driving Force is "
                "not defined — there are no shared intermediates to trade off. "
                + (f"ΔrG′ at 1 mM reference concentrations is <b>{report.mdf:.2f} "
                   "kJ/mol</b>." if report.mdf is not None else "")))
        else:
            verdict = ("<span style='color:#188038'>feasible</span>"
                       if report.mdf_feasible else
                       "<span style='color:#C5221F'>infeasible</span>")
            v.addWidget(_body(
                f"Max-min Driving Force: <b>{report.mdf:.2f} kJ/mol</b> — {verdict}. "
                f"Scored {report.mdf_scored} reaction(s)"
                + (f", {report.mdf_missing} without ΔrG′°" if report.mdf_missing else "")
                + "."))
        if report.dg_per_step:
            rows = "<br>".join(f"&nbsp;&nbsp;<b>{k}</b>: ΔrG′ {val:+.1f} kJ/mol"
                               for k, val in list(report.dg_per_step.items())[:12])
            v.addWidget(_body(rows))
        if report.mdf is not None:
            v.addWidget(_body(
                "Assumed physiology: pH 7.5, ionic strength 0.25 M, pMg 3.0, metabolite "
                "concentrations 1 µM – 10 mM, 298 K.", style.TEXT_MUTED))

        # ---- flux and blockage -------------------------------------------------
        v.addWidget(_section("Flux"))
        if report.blocked_step:
            v.addWidget(_body(
                f"✗ <b>Blocked at {report.blocked_step}</b> — this step cannot carry flux "
                "in any steady state.", "#C5221F"))
            if report.blocking_metabolites:
                v.addWidget(_body("Blocking metabolite(s): <b>"
                                  + ", ".join(report.blocking_metabolites[:6]) + "</b>."))
            if report.recommendation:
                v.addWidget(_body(f"<u>What to do</u><br>{report.recommendation}"))
        elif report.max_flux is None:
            v.addWidget(_body("Flux was not computed for this route.", style.TEXT_MUTED))
        elif report.flux_is_indicative:
            v.addWidget(_body(
                ("✓ This rule-based route <b>can carry flux</b>." if report.max_flux > 1e-9
                 else "✗ This rule-based route <b>carries no flux</b> as written.")
                + " The absolute value is not a validated capacity for rule-derived "
                  "chemistry, so it is not shown."))
        else:
            if report.carbon_yield is not None:
                cy = f" · carbon yield <b>{report.carbon_yield * 100:.1f}%</b> of consumed C"
            elif report.carbon_yield_note:
                cy = (" · carbon yield <b>not computable</b> — "
                      f"{report.carbon_yield_note}")
            else:
                cy = ""
            v.addWidget(_body(
                f"Maximum production: <b>{report.max_flux:.4g}</b> mmol gDW⁻¹ h⁻¹{cy}."))

        # ---- competition -------------------------------------------------------
        v.addWidget(_section("Competition for intermediates"))
        if report.n_competitors:
            v.addWidget(_body(
                f"{report.n_competitors} native reaction(s) compete for this route's "
                "intermediates — these are your knock-down candidates."))
            for met, others in list(report.competitors.items())[:6]:
                v.addWidget(_body(f"&nbsp;&nbsp;<b>{met}</b> is also consumed by: "
                                  + ", ".join(others)))
            if report.linear_yield is not None and report.network_yield is not None:
                v.addWidget(_body(
                    f"Ideal (linear) yield {report.linear_yield:.4g} vs realistic "
                    f"(network) yield {report.network_yield:.4g}."))
        else:
            v.addWidget(_body("✓ No native reaction competes for this route's "
                              "intermediates.", "#188038"))

        scroll.setWidget(inner)
        outer.addWidget(scroll, 1)

        # ---- actions -----------------------------------------------------------
        row = QHBoxLayout()
        if on_fetch_missing is not None and blockage is not None and blockage.blockers:
            btn = QPushButton("Fetch missing reactions…")
            btn.setObjectName("primary")
            btn.setToolTip("Fetch chemistry around the blocking metabolite from KEGG and "
                           "Rhea, then search again for a complete route.")
            btn.clicked.connect(self._fetch)
            row.addWidget(btn)
        row.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        row.addWidget(close)
        outer.addLayout(row)

    def _fetch(self) -> None:
        self.accept()
        if self._on_fetch:
            self._on_fetch()
