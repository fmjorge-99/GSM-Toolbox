"""Mass/charge balance assessment and proton/water balancing.

Shared by the reaction editors and the pathway-design suggestions so the same,
carefully-scoped logic is used everywhere (see the design note below).

Design (grounded in standard GEM curation — Thiele & Palsson 2010; memote;
ModelSEED/BiGG balancing conventions):

* A metabolic (non-boundary) reaction should conserve every element and total
  charge. We only *trust* a verdict when every metabolite has a formula, so we
  first try to back-fill missing formulas from annotation before judging.
* Proton (H⁺) and water (H₂O) balancing is legitimate *only* when the residual is
  confined to H, O and charge — the protonation/hydration bookkeeping curators do
  routinely for a dilute aqueous cytosol at ~pH 7. It couples the reaction to the
  model's already-balanced water/proton pools, so it does not create free mass.
* A residual involving any other element (C, N, P, S, …) signals a genuinely
  wrong/incomplete reaction (a missing substrate or a bad formula) and must NOT be
  auto-patched — we refuse and say so, rather than hiding the error with water.

**Absent is not zero.** cobra's own ``Reaction.check_mass_balance`` reads a metabolite
with no formula as contributing no atoms, so a reaction whose product has no formula is
reported as *losing* that product's mass — ``Fucoxanthin_synthase`` comes back as
``{'C': -42, 'H': -58, 'O': -6}`` when nothing is wrong with the reaction at all; the
product's atoms simply were never counted. Every verdict in this module therefore has
three outcomes, not two: **balanced**, **imbalanced**, and **unverifiable** (a named
participant has no formula). :class:`BalanceVerdict` is the type that carries them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

_FORMULA_ANNOTATION_KEYS = {
    "formula", "chemicalformula", "chemical_formula", "formulae", "seed.formula",
}


def parse_formula(formula: str) -> Dict[str, float]:
    """Element -> count from a chemical formula (e.g. 'C6H12O6' -> {C:6,H:12,O:6})."""
    counts: Dict[str, float] = {}
    for elem, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        if not elem:
            continue
        counts[elem] = counts.get(elem, 0.0) + (float(num) if num else 1.0)
    return counts


def formula_from_annotation(met) -> str:
    """Best-effort formula string from a metabolite's annotation."""
    ann = getattr(met, "annotation", None) or {}
    if isinstance(ann, dict):
        for key, value in ann.items():
            if str(key).lower() in _FORMULA_ANNOTATION_KEYS:
                v = value[0] if isinstance(value, (list, tuple)) else value
                v = str(v).strip()
                if v and v.lower() not in ("nan", "none"):
                    return v
    return ""


def elements_and_charge(met) -> Tuple[Dict[str, float], float, bool]:
    """``(elements, charge, known)`` for a metabolite, back-filling the formula
    from annotation when the attribute is empty.

    ``charge`` is ``None`` when the metabolite carries none. That distinction matters:
    reading an absent charge as 0 makes a reaction look charge-unbalanced whenever some
    of its participants are annotated and others are not — a reaction that is perfectly
    mass-balanced then gets reported as broken, purely because of a gap in the database.
    """
    formula = getattr(met, "formula", None) or formula_from_annotation(met)
    if not formula:
        return {}, None, False
    charge = getattr(met, "charge", None)
    return parse_formula(formula), (None if charge is None else float(charge)), True


def residual_from_terms(terms, tol: float = 1e-6):
    """Net element/charge residual (products − reactants) for ``terms`` =
    iterable of ``(coefficient, metabolite)``. Returns ``(residual, unknown)``
    where ``unknown`` lists metabolites with no usable formula.

    Charge is only reported when *every* participant declares one — a partial charge
    sum is not a balance, it is an artefact of what happens to be annotated.
    """
    net: Dict[str, float] = {}
    unknown = []
    charge_net = 0.0
    charge_known = True
    for coeff, met in terms:
        els, charge, known = elements_and_charge(met)
        if not known:
            unknown.append(met)
            continue
        for e, n in els.items():
            net[e] = net.get(e, 0.0) + coeff * n
        if charge is None:
            charge_known = False
        else:
            charge_net += coeff * charge
    if charge_known:
        net["charge"] = charge_net
    net = {k: v for k, v in net.items() if abs(v) > tol}
    return net, unknown


#: The three possible outcomes of a balance check. ``UNVERIFIABLE`` is a statement about
#: our *data*, not about the reaction: it must never be shown as, or counted with, a
#: genuine imbalance.
BALANCED = "balanced"
IMBALANCED = "imbalanced"
UNVERIFIABLE = "unverifiable"


def metabolite_label(met) -> str:
    """``'name (id)'`` for a metabolite, or just the id when there is no distinct name."""
    mid = str(getattr(met, "id", met) or "")
    name = str(getattr(met, "name", "") or "").strip()
    return f"{name} ({mid})" if name and name != mid else mid


@dataclass(frozen=True)
class BalanceVerdict:
    """The outcome of checking one reaction's mass/charge balance.

    ``status`` is one of :data:`BALANCED`, :data:`IMBALANCED`, :data:`UNVERIFIABLE`.
    ``residual`` is only meaningful for :data:`IMBALANCED`; for :data:`UNVERIFIABLE` it is
    empty, because a residual computed over some of the participants is not a residual —
    it is the mass of whatever we failed to count.
    """

    status: str = BALANCED
    residual: Dict[str, float] = field(default_factory=dict)
    #: ids of the participants that carry no usable formula (empty unless UNVERIFIABLE)
    missing_formula: Tuple[str, ...] = ()
    #: the same participants, as ``'name (id)'`` for display
    missing_formula_labels: Tuple[str, ...] = ()

    @property
    def checkable(self) -> bool:
        return self.status != UNVERIFIABLE

    @property
    def balanced(self) -> bool:
        """True only when the reaction was checked *and* found balanced.

        Deliberately False for an unverifiable reaction: callers that gate on
        ``balanced`` are asking "may I trust this?", and the honest answer is no. Callers
        that must not *penalise* a data gap test ``checkable`` first — see
        :attr:`verifiably_imbalanced`.
        """
        return self.status == BALANCED

    @property
    def verifiably_imbalanced(self) -> bool:
        """True only for a real imbalance, with every participant's formula known."""
        return self.status == IMBALANCED

    def as_tuple(self) -> Tuple[bool, Dict[str, float], bool]:
        """``(balanced, residual, checkable)`` — the older three-value contract.

        ``balanced`` is reported as True for an unverifiable reaction here, because the
        callers of this form all guard with ``checkable`` and read the pair as "nothing
        found against it".
        """
        return (self.status != IMBALANCED, dict(self.residual), self.checkable)

    def sentence(self, *, max_named: int = 3) -> str:
        """A one-line verdict fit to show a user."""
        if self.status == UNVERIFIABLE:
            return f"cannot be checked — {self.missing_formula_sentence(max_named=max_named)}"
        if self.status == IMBALANCED:
            return f"mass/charge unbalanced — {format_residual(self.residual)}"
        return "mass and charge balanced"

    def missing_formula_sentence(self, *, max_named: int = 3) -> str:
        """``'Fucoxanthin (Fucoxanthin_c) has no formula'`` — names the actual culprit."""
        labels = list(self.missing_formula_labels)
        if not labels:
            return ""
        shown = ", ".join(labels[:max_named])
        extra = len(labels) - max_named
        if extra > 0:
            shown += f" and {extra} other(s)"
        verb = "has" if len(labels) == 1 else "have"
        return f"{shown} {verb} no formula"


def assess_terms(terms: Iterable[Tuple[float, object]]) -> BalanceVerdict:
    """Balance verdict for ``terms`` = iterable of ``(coefficient, metabolite)``."""
    net, unknown = residual_from_terms(terms)
    if unknown:
        return BalanceVerdict(
            status=UNVERIFIABLE,
            missing_formula=tuple(str(getattr(m, "id", m)) for m in unknown),
            missing_formula_labels=tuple(metabolite_label(m) for m in unknown),
        )
    return BalanceVerdict(status=(IMBALANCED if net else BALANCED), residual=net)


def assess_reaction(reaction) -> BalanceVerdict:
    """Balance verdict for a cobra Reaction.

    Prefer this over :func:`reaction_residual`: it distinguishes "this reaction does not
    balance" from "we have no formula for participant X", which the residual alone
    cannot express.
    """
    return assess_terms((c, m) for m, c in reaction.metabolites.items())


def reaction_residual(reaction):
    """Element/charge residual of a cobra Reaction, or ``None`` if any metabolite
    has no usable formula (so the reaction can't be verified).

    ``None`` (unverifiable) and ``{}`` (balanced) are different answers — a caller
    writing ``if not residual`` silently merges them. Use :func:`assess_reaction` when
    the difference matters, which is nearly always.
    """
    verdict = assess_reaction(reaction)
    return None if not verdict.checkable else verdict.residual


def format_residual(residual: dict) -> str:
    return ", ".join(f"{k}{v:+.3g}" for k, v in sorted(residual.items()))


def missing_formula_metabolites(reaction) -> List:
    """The reaction's participants that carry no usable formula (possibly empty)."""
    _net, unknown = residual_from_terms((c, m) for m, c in reaction.metabolites.items())
    return unknown


def water_proton_fix(residual: dict) -> Optional[Tuple[float, float]]:
    """``(n_water, n_proton)`` closing ``residual`` with only H₂O/H⁺, or ``None``
    if it can't be (other elements present, or an inconsistent H/O/charge residual).
    Positive counts belong on the product side, negative on the reactant side."""
    if set(residual) - {"H", "O", "charge"}:
        return None
    n_w = -residual.get("O", 0.0)
    n_h = -residual.get("charge", 0.0)
    if abs(residual.get("H", 0.0) + 2 * n_w + n_h) > 1e-6:
        return None
    if abs(n_w) < 1e-9 and abs(n_h) < 1e-9:
        return None
    return n_w, n_h


def dominant_compartment(reaction) -> str:
    from collections import Counter
    cnt = Counter(m.compartment or "c" for m in reaction.metabolites)
    return cnt.most_common(1)[0][0] if cnt else "c"


def plan_water_proton_changes(reaction) -> Optional[List[dict]]:
    """Plan (don't apply) the H₂O/H⁺ additions that would balance ``reaction``.

    Returns a list of ``{metabolite_id, delta, side}`` changes, or ``None`` if the
    reaction is already balanced, unverifiable, or not fixable with water/protons.
    ``delta`` is the signed stoichiometric coefficient to add (products +, reactants −).
    """
    verdict = assess_reaction(reaction)
    if not verdict.verifiably_imbalanced:
        return None            # already balanced, or a formula is missing: nothing to fix
    fix = water_proton_fix(verdict.residual)
    if fix is None:
        return None
    n_w, n_h = fix
    comp = dominant_compartment(reaction)
    model = reaction.model
    changes = []
    for base, amount in (("h2o", n_w), ("h", n_h)):
        if abs(amount) < 1e-9:
            continue
        mid = f"{base}_{comp}"
        if model is None or not model.metabolites.has_id(mid):
            return None    # need real water/proton metabolites in the model
        changes.append({
            "metabolite_id": mid,
            "delta": float(amount),
            "side": "product" if amount > 0 else "reactant",
        })
    return changes or None


def apply_changes(reaction, changes: List[dict]) -> None:
    """Apply planned changes (from :func:`plan_water_proton_changes`) to a cobra
    Reaction in place, merging with any existing coefficient."""
    model = reaction.model
    add = {}
    for ch in changes:
        met = model.metabolites.get_by_id(ch["metabolite_id"])
        add[met] = add.get(met, 0.0) + ch["delta"]
    reaction.add_metabolites(add)
