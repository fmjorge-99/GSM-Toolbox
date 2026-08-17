"""Explore heterologous chemistry *upstream* of a designed route's entry point.

A finished design ends where it meets the host: some native compound the route draws on
as its precursor. That hand-over point is the least examined part of a design, and it is
often where the yield is actually decided —

* the precursor may be **idle** (present in the model but carrying no flux), so the route
  cannot run at all until something supplies it. Lactaldehyde in Synechocystis is the
  standing example: a 1,2-propanediol design that starts there is not buildable until a
  methylglyoxal synthase is added upstream;
* or the precursor carries flux, but only barely, and a heterologous step feeding it
  would lift the whole route.

Either way the question is the same: **what else could make this compound?** This module
searches the loaded databases for reactions that produce a route's entry metabolites, and
ranks them by whether their own substrates are actually available in the host.

It is strictly additive — it reads a finished route and proposes extra reactions. Nothing
here changes how a route is found.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import cobra

from . import flux_context as fx


@dataclass
class UpstreamCandidate:
    """One heterologous reaction that could feed a route's entry metabolite."""

    reaction_id: str
    name: str = ""
    equation: str = ""
    ec_numbers: List[str] = field(default_factory=list)
    #: Substrates it needs, split by whether the host can actually supply them.
    available_substrates: List[str] = field(default_factory=list)
    missing_substrates: List[str] = field(default_factory=list)
    idle_substrates: List[str] = field(default_factory=list)
    balanced: bool = True
    balance_checkable: bool = True

    @property
    def ready(self) -> bool:
        """True when every substrate is already produced by the host at growth."""
        return not self.missing_substrates and not self.idle_substrates

    @property
    def score(self) -> tuple:
        """Sort key — best first. Reactions the host could run today come first."""
        return (len(self.missing_substrates), len(self.idle_substrates),
                0 if self.balanced else 1, self.reaction_id)

    def verdict(self) -> str:
        if self.missing_substrates:
            return (f"needs {len(self.missing_substrates)} compound(s) the host does not "
                    "have: " + ", ".join(self.missing_substrates[:3]))
        if self.idle_substrates:
            return ("all substrates exist but "
                    + ", ".join(self.idle_substrates[:3]) + " carry no flux — this step "
                    "would need its own upstream supply")
        return "every substrate is actively produced by the host"


@dataclass
class UpstreamReport:
    """Entry points of a route and what could feed each of them."""

    entry_metabolites: Dict[str, str] = field(default_factory=dict)   # id -> name
    idle_entries: Dict[str, str] = field(default_factory=dict)        # id -> name
    candidates: Dict[str, List[UpstreamCandidate]] = field(default_factory=dict)
    context_summary: str = ""

    @property
    def total_candidates(self) -> int:
        return sum(len(v) for v in self.candidates.values())

    def headline(self) -> str:
        if not self.entry_metabolites:
            return "This route draws nothing from the host — there is no entry point."
        n = len(self.entry_metabolites)
        if self.idle_entries:
            names = ", ".join(list(self.idle_entries.values())[:3])
            return (f"This route enters host metabolism at {n} compound(s). "
                    f"<b>{len(self.idle_entries)} carry no flux</b> ({names}) — the route "
                    "cannot run until something upstream supplies them.")
        return (f"This route enters host metabolism at {n} compound(s), all of which the "
                "host actively produces.")


def entry_metabolites(universal: cobra.Model, reaction_ids: List[str],
                      host: cobra.Model) -> Dict[str, str]:
    """Compounds the route consumes but does not make itself, and that the host has.

    These are the hand-over points where the design meets native metabolism.
    """
    consumed: Set[str] = set()
    produced: Set[str] = set()
    for rid in reaction_ids:
        if not universal.reactions.has_id(rid):
            continue
        for met, coeff in universal.reactions.get_by_id(rid).metabolites.items():
            (consumed if coeff < 0 else produced).add(met.id)
    from .pathway_design import _is_currency_met

    out: Dict[str, str] = {}
    for mid in sorted(consumed - produced):
        if not universal.metabolites.has_id(mid):
            continue
        met = universal.metabolites.get_by_id(mid)
        if _is_currency_met(met):
            continue                      # ATP/NADH are not the interesting hand-over
        if host.metabolites.has_id(mid):
            out[mid] = met.name or mid
    return out


def explore_upstream(host: cobra.Model, universal: cobra.Model,
                     reaction_ids: List[str], *,
                     limit_per_metabolite: int = 25,
                     context: Optional[fx.FluxContext] = None) -> UpstreamReport:
    """Find database reactions that could feed this route's entry metabolites."""
    report = UpstreamReport()
    ctx = context or fx.growth_flux_context(host)
    report.context_summary = ctx.summary()

    report.entry_metabolites = entry_metabolites(universal, reaction_ids, host)
    if ctx.status == "optimal":
        report.idle_entries = {mid: name
                               for mid, name in report.entry_metabolites.items()
                               if mid not in ctx.produced}

    route = set(reaction_ids)
    host_ids = {m.id for m in host.metabolites}
    from .pathway_design import _is_currency_met, reaction_balance

    for mid in report.entry_metabolites:
        if not universal.metabolites.has_id(mid):
            continue
        met = universal.metabolites.get_by_id(mid)
        found: List[UpstreamCandidate] = []
        for rxn in met.reactions:
            if rxn.id in route or rxn.boundary:
                continue
            coeff = rxn.metabolites[met]
            # Producing it as written, or in reverse if the reaction is reversible.
            makes = (coeff > 0 and rxn.upper_bound > 1e-9) or \
                    (coeff < 0 and rxn.lower_bound < -1e-9)
            if not makes:
                continue
            forward = coeff > 0
            cand = UpstreamCandidate(
                reaction_id=rxn.id,
                name=rxn.name or "",
                equation=rxn.build_reaction_string(use_metabolite_names=True))
            annotation = getattr(rxn, "annotation", {}) or {}
            ec = annotation.get("ec-code") or annotation.get("EC Number")
            if ec:
                cand.ec_numbers = ec if isinstance(ec, list) else [str(ec)]
            for other, c in rxn.metabolites.items():
                if other.id == mid or _is_currency_met(other):
                    continue
                # Substrates are the ones consumed in the direction that makes our target.
                is_substrate = (c < 0) if forward else (c > 0)
                if not is_substrate:
                    continue
                if other.id not in host_ids:
                    cand.missing_substrates.append(other.name or other.id)
                elif ctx.status == "optimal" and other.id not in ctx.produced:
                    cand.idle_substrates.append(other.name or other.id)
                else:
                    cand.available_substrates.append(other.name or other.id)
            balanced, _residual, checkable = reaction_balance(rxn)
            cand.balanced, cand.balance_checkable = balanced, checkable
            found.append(cand)
        found.sort(key=lambda c: c.score)
        if found:
            report.candidates[mid] = found[:limit_per_metabolite]
    return report
