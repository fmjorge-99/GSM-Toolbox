"""A single, honest verdict on whether a designed pathway is actually buildable.

The results table makes it far too easy to rank routes by predicted flux, which is the
one number that most often points at the wrong route: across the compounds studied, the
highest-flux route was repeatedly mislabelled chemistry, thermodynamically impossible, or
smothered by competing reactions. This module gathers everything the toolbox knows about
a route — chemistry, balance, thermodynamics, competition, flux, blockage — into one
sentence the user reads first, and one detailed report they can open.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

# Verdict levels, worst first when summarising.
BLOCKED = "blocked"            # cannot carry flux at all
IMPLAUSIBLE = "implausible"    # chemistry looks wrong (isomer/backbone problem)
INFEASIBLE = "infeasible"      # thermodynamically impossible
CAUTION = "caution"            # runs, but with a caveat worth reading
GOOD = "good"                  # nothing known against it

_LABEL = {
    BLOCKED: "Blocked",
    IMPLAUSIBLE: "Chemistry questionable",
    INFEASIBLE: "Thermodynamically infeasible",
    CAUTION: "Usable with caveats",
    GOOD: "Looks buildable",
}
_COLOUR = {
    BLOCKED: "#C5221F", IMPLAUSIBLE: "#C5221F", INFEASIBLE: "#E8710A",
    CAUTION: "#E8710A", GOOD: "#188038",
}


@dataclass
class FeasibilityReport:
    """Everything known about one designed route."""

    target: str = ""
    reaction_ids: List[str] = field(default_factory=list)
    verdict: str = GOOD

    # flux
    max_flux: Optional[float] = None
    carbon_yield: Optional[float] = None
    #: why there is no carbon yield, when there is none — "not computable because the
    #: product has no formula" is a fixable database gap, not a failed design
    carbon_yield_note: str = ""
    flux_is_indicative: bool = False

    # chemistry
    balanced: bool = True
    unverified_steps: List[str] = field(default_factory=list)
    #: step id -> the participants whose missing formula blocked its balance check
    unverified_reasons: Dict[str, List[str]] = field(default_factory=dict)
    isomer_checked: bool = False
    isomer_coverage: str = ""
    isomer_warnings: List[str] = field(default_factory=list)
    #: what the automatic chemistry back-fill recovered before the checks ran
    backfill_note: str = ""

    # thermodynamics
    mdf: Optional[float] = None
    mdf_feasible: Optional[bool] = None
    mdf_single_reaction: bool = False
    mdf_scored: int = 0
    mdf_missing: int = 0
    mdf_note: str = ""
    dg_per_step: Dict[str, float] = field(default_factory=dict)

    # network
    n_competitors: int = 0
    competitors: Dict[str, List[str]] = field(default_factory=dict)
    linear_yield: Optional[float] = None
    network_yield: Optional[float] = None

    # blockage
    blocked_step: str = ""
    blocking_metabolites: List[str] = field(default_factory=list)
    recommendation: str = ""

    @property
    def label(self) -> str:
        return _LABEL.get(self.verdict, self.verdict)

    @property
    def colour(self) -> str:
        return _COLOUR.get(self.verdict, "#5f6368")

    def sentence(self) -> str:
        """The one-line verdict shown beside every result (VI.10)."""
        bits: List[str] = []

        # chemistry
        if self.isomer_warnings:
            bits.append("chemistry questionable")
        elif self.isomer_checked:
            bits.append("isomer-checked ✓")
        else:
            bits.append("isomers not checkable")

        # balance. "Unverifiable" is a third outcome, not a shade of "unbalanced": the
        # participant's atoms were never counted, so nothing was found against the step.
        if not self.balanced:
            bits.append("unbalanced step(s)")
        elif self.unverified_steps:
            bits.append(f"balanced ✓ ({len(self.unverified_steps)} unverifiable — "
                        "missing formula)")
        else:
            bits.append("balanced ✓")

        # thermodynamics
        if self.mdf_feasible is True:
            bits.append("thermodynamically feasible ✓" if not self.mdf_single_reaction
                        else "favourable ✓ (single reaction)")
        elif self.mdf_feasible is False:
            bits.append("thermodynamically infeasible ✗")
        else:
            bits.append("thermodynamics not available")

        # flux
        if self.blocked_step:
            bits.append(f"blocked at {self.blocked_step}")
        elif self.max_flux is None:
            bits.append("flux unknown")
        elif self.flux_is_indicative:
            bits.append("carries flux" if self.max_flux > 1e-9 else "carries no flux")
        elif self.max_flux > 1e-9:
            if self.carbon_yield is not None and self.carbon_yield == self.carbon_yield:
                cy = f", {self.carbon_yield * 100:.0f}% C-yield"
            elif self.carbon_yield_note:
                cy = f", C-yield not computable ({self.carbon_yield_note})"
            else:
                cy = ""
            bits.append(f"flux {self.max_flux:.3g}{cy}")
        else:
            bits.append("no flux")

        # competition
        if self.n_competitors:
            bits.append(f"{self.n_competitors} competing reaction(s)")
        else:
            bits.append("no competitors")

        return f"{self.label} — " + "; ".join(bits) + "."


def _worst(*verdicts: str) -> str:
    order = [BLOCKED, IMPLAUSIBLE, INFEASIBLE, CAUTION, GOOD]
    for v in order:
        if v in verdicts:
            return v
    return GOOD


def assess(result, *, diagnosis=None, branching=None, mdf=None) -> FeasibilityReport:
    """Build the report from a :class:`PathwayResult` plus any analyses already run.

    Every input is optional: the report degrades gracefully and says which parts are
    unknown rather than implying a clean bill of health.
    """
    r = FeasibilityReport(
        target=getattr(result, "target", ""),
        reaction_ids=list(getattr(result, "reaction_ids", []) or []),
        balanced=bool(getattr(result, "balanced", True)),
        unverified_steps=list(getattr(result, "unverified_steps", []) or []),
        unverified_reasons=dict(getattr(result, "unverified_reasons", None) or {}),
        isomer_checked=bool(getattr(result, "isomer_checked", False)),
        isomer_coverage=getattr(result, "isomer_coverage", "") or "",
        isomer_warnings=list(getattr(result, "isomer_warnings", []) or []),
        flux_is_indicative=bool(getattr(result, "flux_is_indicative", False)),
    )
    f = getattr(result, "production_flux", float("nan"))
    r.max_flux = None if f != f else float(f)
    cy = getattr(result, "carbon_yield", float("nan"))
    r.carbon_yield = None if cy != cy else float(cy)
    r.carbon_yield_note = "" if r.carbon_yield is not None else \
        (getattr(result, "carbon_yield_note", "") or "")

    if diagnosis is not None:
        r.blocked_step = getattr(diagnosis, "bottleneck", "") or ""
        r.blocking_metabolites = list(getattr(diagnosis, "blocking_metabolites", []) or [])
        r.recommendation = getattr(diagnosis, "recommendation", "") or ""
    if branching is not None:
        r.n_competitors = int(getattr(branching, "n_effective_branches", 0) or 0)
        r.competitors = {k: list(v)[:6] for k, v in
                         (getattr(branching, "per_intermediate", {}) or {}).items()}
        ly = getattr(branching, "linear_yield", float("nan"))
        ny = getattr(branching, "network_yield", float("nan"))
        r.linear_yield = None if ly != ly else float(ly)
        r.network_yield = None if ny != ny else float(ny)
    if mdf is not None:
        m = getattr(mdf, "mdf", float("nan"))
        r.mdf = None if m != m else float(m)
        r.mdf_feasible = bool(getattr(mdf, "feasible", False)) if r.mdf is not None else None
        r.mdf_single_reaction = bool(getattr(mdf, "single_reaction", False))
        r.mdf_scored = len(getattr(mdf, "dg0_prime", {}) or {})
        r.mdf_missing = len(getattr(mdf, "missing", []) or [])
        r.mdf_note = getattr(mdf, "note", "") or ""
        r.dg_per_step = dict(getattr(mdf, "dg_prime", {}) or {})

    # ---- overall verdict --------------------------------------------------------
    verdicts = []
    if r.isomer_warnings:
        verdicts.append(IMPLAUSIBLE)
    if r.mdf_feasible is False:
        verdicts.append(INFEASIBLE)
    if r.blocked_step or (r.max_flux is not None and r.max_flux <= 1e-9
                          and not r.flux_is_indicative):
        verdicts.append(BLOCKED)
    if not r.balanced or r.n_competitors:
        verdicts.append(CAUTION)
    r.verdict = _worst(*verdicts) if verdicts else GOOD
    return r
