"""A managed library of regulatory rule sets, mirroring how reaction databases work.

Regulation is organism-specific in a way that reaction chemistry is not. A rule saying
"NdhR de-represses the CCM below 1 mM inorganic carbon" is a statement about
*Synechocystis*; applied to an *E. coli* model it is not merely unhelpful, it is wrong
while looking authoritative. Shipping any rule set as a built-in default would therefore
attach cyanobacterial physiology to whatever model happened to be open.

So rule sets are **files the user loads**, not code:

* nothing is active until the user chooses a set — with no choice, ``active()`` returns an
  empty rule set, and an empty rule set reproduces the unregulated model exactly;
* a loaded file is copied into the library under ``~/.gsm_toolbox/regulation`` and stays
  available across sessions, the same contract as a downloaded reaction database;
* every file is validated before it is trusted, and the problems are reported rather than
  swallowed — a malformed rule that silently did nothing would be worse than one that
  refuses to load.

The *Synechocystis* set developed in this project is still available, but as an **example
to import**, never as a default (:func:`examples`).
"""

from __future__ import annotations

import json
import os
import shutil
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import cache
from . import regulation as reg


def library_dir() -> str:
    path = os.path.join(cache.base_dir(), "regulation")
    os.makedirs(path, exist_ok=True)
    return path


@dataclass
class RulesetInfo:
    """What the manager shows about one stored rule set."""

    name: str
    path: str
    n_rules: int = 0
    n_enabled: int = 0
    organism: str = ""
    description: str = ""
    source: str = ""
    imported_at: float = 0.0
    problems: List[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.problems

    def imported_text(self) -> str:
        if not self.imported_at:
            return "—"
        return time.strftime("%Y-%m-%d %H:%M", time.localtime(self.imported_at))

    def weakest_confidence(self) -> str:
        """Worst confidence among the enabled rules, or ``measured`` if none."""
        try:
            ruleset = reg.load(self.path)
        except Exception:  # noqa: BLE001 — a broken file is reported via `problems`
            return reg.MEASURED
        levels = [r.confidence for r in ruleset.enabled()]
        if not levels:
            return reg.MEASURED
        order = {reg.MEASURED: 0, reg.INFERRED: 1, reg.ASSUMED: 2}
        return max(levels, key=lambda c: order.get(c, 9))


# --------------------------------------------------------------------------------------
# Validation — a user-supplied file is untrusted input
# --------------------------------------------------------------------------------------
_REQUIRED = ("name", "sensor", "effect")
_EFFECTS = {reg.GATE, reg.SCALE, reg.ENZYME_COST, reg.BUDGET, reg.BIOMASS, reg.PARAMETER}
_CONFIDENCES = {reg.MEASURED, reg.INFERRED, reg.ASSUMED}


def validate(payload: dict) -> List[str]:
    """Return every problem found in a rule-set payload; empty means usable.

    Reports *all* problems rather than the first, so a user fixing a hand-edited file can
    see the whole list instead of discovering them one reload at a time.
    """
    problems: List[str] = []
    if not isinstance(payload, dict):
        return ["The file does not contain a rule-set object."]
    rules = payload.get("rules")
    if not isinstance(rules, list):
        return ["The file has no 'rules' list."]

    seen = set()
    for index, item in enumerate(rules, start=1):
        where = f"rule {index}"
        if not isinstance(item, dict):
            problems.append(f"{where}: not an object")
            continue
        name = item.get("name")
        if name:
            where = f"'{name}'"
            if name in seen:
                problems.append(f"{where}: duplicate rule name")
            seen.add(name)
        for key in _REQUIRED:
            if not item.get(key):
                problems.append(f"{where}: missing '{key}'")
        effect = item.get("effect")
        if effect and effect not in _EFFECTS:
            problems.append(f"{where}: unknown effect '{effect}' "
                            f"(expected one of {', '.join(sorted(_EFFECTS))})")
        confidence = item.get("confidence", reg.ASSUMED)
        if confidence not in _CONFIDENCES:
            problems.append(f"{where}: unknown confidence '{confidence}'")
        spec = item.get("response") or {}
        kind = spec.get("kind")
        if kind not in reg.RESPONSES:
            problems.append(f"{where}: unknown response '{kind}' "
                            f"(expected one of {', '.join(sorted(reg.RESPONSES))})")
        else:
            try:
                reg.RESPONSES[kind](**{k: v for k, v in spec.items() if k != "kind"})
            except TypeError as exc:
                problems.append(f"{where}: response parameters do not fit '{kind}' ({exc})")
        if effect in (reg.GATE, reg.SCALE, reg.ENZYME_COST, reg.BIOMASS, reg.PARAMETER) \
                and not item.get("targets"):
            problems.append(f"{where}: '{effect}' needs at least one target")
    return problems


def inspect(path: str) -> RulesetInfo:
    """Read one rule-set file into a summary, without trusting its contents."""
    info = RulesetInfo(name=os.path.splitext(os.path.basename(path))[0], path=path)
    try:
        with open(path, encoding="utf-8") as fh:
            payload = json.load(fh)
    except Exception as exc:  # noqa: BLE001 — unreadable is a reportable state
        info.problems.append(f"cannot be read: {exc}")
        return info

    info.problems = validate(payload)
    rules = payload.get("rules") if isinstance(payload, dict) else []
    if isinstance(rules, list):
        info.n_rules = len(rules)
        info.n_enabled = sum(1 for r in rules
                             if isinstance(r, dict) and r.get("enabled", True))
    if isinstance(payload, dict):
        info.name = payload.get("name") or info.name
        info.organism = payload.get("organism", "")
        info.description = payload.get("description", "")
        info.source = payload.get("source", "")
    try:
        info.imported_at = os.path.getmtime(path)
    except OSError:
        pass
    return info


# --------------------------------------------------------------------------------------
# The library
# --------------------------------------------------------------------------------------
def scan_library() -> List[RulesetInfo]:
    """Every rule set stored locally, newest first."""
    out = []
    for entry in sorted(os.listdir(library_dir())):
        if entry.lower().endswith(".json"):
            out.append(inspect(os.path.join(library_dir(), entry)))
    out.sort(key=lambda i: i.imported_at, reverse=True)
    return out


def _unique_path(stem: str) -> str:
    safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in stem).strip() or "rules"
    path = os.path.join(library_dir(), f"{safe}.json")
    n = 2
    while os.path.exists(path):
        path = os.path.join(library_dir(), f"{safe} ({n}).json")
        n += 1
    return path


def import_file(source_path: str) -> RulesetInfo:
    """Copy a rule-set file into the library after validating it.

    An invalid file is *not* imported. Storing something the app will later refuse to use
    would leave the user with a library entry that fails at the point of simulation, far
    from the action that caused it.
    """
    info = inspect(source_path)
    if not info.valid:
        return info
    stem = info.name or os.path.splitext(os.path.basename(source_path))[0]
    destination = _unique_path(stem)
    shutil.copyfile(source_path, destination)
    return inspect(destination)


def store(ruleset: reg.RuleSet, name: str = "", path: str = "") -> RulesetInfo:
    """Write a rule set into the library (or back to ``path`` if it is already there)."""
    target = path or _unique_path(name or ruleset.name or "rules")
    reg.save(ruleset, target)
    return inspect(target)


def remove(path: str) -> bool:
    """Delete a stored rule set. Clears the active selection if it was the active one."""
    try:
        os.remove(path)
    except OSError:
        return False
    if os.path.abspath(active_path() or "") == os.path.abspath(path):
        set_active("")
    return True


# --------------------------------------------------------------------------------------
# Examples — offered, never applied by default
# --------------------------------------------------------------------------------------
def examples() -> List[Tuple[str, str]]:
    """``(label, path)`` for the rule sets shipped as importable examples.

    These are *offered*, not active. The distinction matters: the Synechocystis set
    encodes cyanobacterial physiology, and silently applying it to another organism's
    model would produce confident, wrong answers.
    """
    from .. import resources

    base = os.path.join(os.path.dirname(os.path.abspath(resources.__file__)),
                        "regulation")
    out = []
    if os.path.isdir(base):
        for entry in sorted(os.listdir(base)):
            if not entry.lower().endswith(".json"):
                continue
            info = inspect(os.path.join(base, entry))
            label = info.name or entry
            if info.organism:
                label = f"{label} — {info.organism}"
            out.append((label, os.path.join(base, entry)))
    return out


# --------------------------------------------------------------------------------------
# Which rule set is in force
# --------------------------------------------------------------------------------------
def active_path() -> str:
    from . import preferences as prefs

    return (prefs.get(prefs.REGULATION_RULESET) or "").strip()


def set_active(path: str) -> None:
    from . import preferences as prefs

    prefs.set(prefs.REGULATION_RULESET, path or "")


def active() -> Tuple[reg.RuleSet, str]:
    """The rule set in force, and where it came from.

    Returns an **empty** rule set when the user has not chosen one — deliberately, so a
    fresh install simulates the plain model rather than somebody else's organism.
    """
    path = active_path()
    if not path or not os.path.exists(path):
        return reg.RuleSet(name=""), ""
    try:
        return reg.load(path), path
    except Exception:  # noqa: BLE001 — a broken active file must not block the app
        return reg.RuleSet(name=""), ""


# --------------------------------------------------------------------------------------
# Model fit — does this rule set actually address the loaded model?
# --------------------------------------------------------------------------------------
@dataclass
class FitReport:
    """How much of a rule set the loaded model can actually be affected by."""

    matched: Dict[str, List[str]] = field(default_factory=dict)   # rule → targets found
    missing: Dict[str, List[str]] = field(default_factory=dict)   # rule → targets absent
    inapplicable: List[str] = field(default_factory=list)         # rules with no target hit

    @property
    def n_targets(self) -> int:
        return sum(len(v) for v in self.matched.values()) + \
               sum(len(v) for v in self.missing.values())

    @property
    def n_matched(self) -> int:
        return sum(len(v) for v in self.matched.values())

    def summary(self) -> str:
        if not self.n_targets:
            return "This rule set names no reaction targets."
        pct = 100.0 * self.n_matched / self.n_targets
        text = (f"{self.n_matched} of {self.n_targets} rule targets exist in this model "
                f"({pct:.0f}%).")
        if self.inapplicable:
            text += (f" {len(self.inapplicable)} rule(s) match nothing and will have no "
                     f"effect: {', '.join(self.inapplicable[:4])}"
                     + (" …" if len(self.inapplicable) > 4 else "") + ".")
        return text


#: Effects that are recorded but never applied to the model.
#:
#: A ``parameter`` rule multiplies a named quantity — the biomass-specific absorption
#: cross-section, for instance — and stores the result in ``RegulationResult.parameters``.
#: Nothing reads that dictionary. The rule therefore appears in every report as having
#: fired while changing no flux and no bound, which is the most misleading state a rule
#: can be in: visible, plausible, and inert. Until a consumer exists, say so.
INERT_EFFECTS = {reg.PARAMETER}


@dataclass
class Effectiveness:
    """Whether each rule can actually change this model's behaviour.

    A rule fails to bite for three quite different reasons, and telling them apart is the
    difference between "the biology does not apply here" and "the file has a typo":

    * **dead** — none of its targets exist in the model, usually an id-convention
      mismatch rather than an intended omission;
    * **inert** — the effect type is not wired into the simulation at all;
    * **slack** — the targets exist but the rule only *raises* a capacity, so it changes
      nothing unless that capacity was limiting.
    """

    dead: Dict[str, List[str]] = field(default_factory=dict)
    partial: Dict[str, List[str]] = field(default_factory=dict)
    inert: List[str] = field(default_factory=list)
    raises_only: List[str] = field(default_factory=list)
    effective: List[str] = field(default_factory=list)

    def warnings(self) -> List[str]:
        out = []
        for name, missing in self.dead.items():
            if missing in (["prot_pool"], ["ER_pool_TG_"]):
                out.append(
                    f"'{name}' reallocates enzyme, but this model is not "
                    f"enzyme-constrained (no {missing[0]}) — it will never do anything.")
            else:
                out.append(f"'{name}' targets nothing in this model "
                           f"({', '.join(missing[:4])}) — it will never do anything.")
        for name, missing in self.partial.items():
            out.append(f"'{name}' is missing {len(missing)} of its targets "
                       f"({', '.join(missing[:3])}).")
        for name in self.inert:
            out.append(f"'{name}' uses an effect the simulation does not apply yet, so "
                       f"it is reported as firing while changing nothing.")
        for name in self.raises_only:
            out.append(f"'{name}' only raises a capacity; it changes the result only if "
                       f"that capacity was limiting.")
        return out

    def summary(self) -> str:
        total = (len(self.effective) + len(self.dead) + len(self.inert)
                 + len(self.raises_only))
        if not total:
            return "No rules to check."
        parts = [f"{len(self.effective)} of {total} rule(s) can change this model"]
        if self.dead:
            parts.append(f"{len(self.dead)} target nothing")
        if self.inert:
            parts.append(f"{len(self.inert)} use an unapplied effect")
        if self.raises_only:
            parts.append(f"{len(self.raises_only)} only raise a capacity")
        return "; ".join(parts) + "."


def effectiveness(ruleset: reg.RuleSet, model) -> Effectiveness:
    """Audit a rule set against a model *before* trusting a simulation that used it.

    Written after a study in which nine rules were reported as firing and exactly two
    changed the answer. Everything else was a dead target, an unapplied effect, or a
    ceiling raised above a flux that never reached it — none of which was visible in the
    output.
    """
    report = Effectiveness()
    # Both protein-allocation effects need an enzyme-constrained model. Applied to a plain
    # stoichiometric reconstruction they run without error and change nothing, which is
    # exactly the failure this audit exists to catch.
    has_pool = model.metabolites.has_id("prot_pool")
    has_budget = model.reactions.has_id("ER_pool_TG_")

    for rule in ruleset.enabled():
        if rule.effect in INERT_EFFECTS:
            report.inert.append(rule.name)
            continue
        if rule.effect == reg.BUDGET:
            # Acts on the whole protein pool; it needs no named target.
            if has_budget:
                report.effective.append(rule.name)
            else:
                report.dead[rule.name] = ["ER_pool_TG_"]
            continue
        if rule.effect == reg.ENZYME_COST and not has_pool:
            report.dead[rule.name] = ["prot_pool"]
            continue

        found, absent = [], []
        for target in rule.targets:
            if rule.effect == reg.BIOMASS:
                exists = model.metabolites.has_id(target)
            else:
                exists = bool(reg._targets_in(model, [target]))
            (found if exists else absent).append(target)

        if not found:
            report.dead[rule.name] = absent or list(rule.targets)
            continue
        if absent:
            report.partial[rule.name] = absent
        # A SCALE rule with magnitude > 1 lifts a bound; it only matters when that bound
        # was binding, which is rarely true in a model limited by something else.
        if rule.effect == reg.SCALE and rule.magnitude > 1.0:
            report.raises_only.append(rule.name)
        else:
            report.effective.append(rule.name)
    return report


def fit(ruleset: reg.RuleSet, model) -> FitReport:
    """Check a rule set against a model *before* it is used to draw conclusions.

    A rule whose targets are absent does not fail — it quietly does nothing, and a scan
    then shows "no regulatory transition" for a reason that has nothing to do with
    biology. Surfacing the mismatch is the difference between a null result and a
    misconfiguration.
    """
    report = FitReport()
    for rule in ruleset.rules:
        if rule.effect in (reg.BUDGET, reg.PARAMETER):
            continue                    # these do not name model reactions
        found, absent = [], []
        for target in rule.targets:
            if rule.effect == reg.BIOMASS:
                exists = model.metabolites.has_id(target)
            else:
                exists = bool(reg._targets_in(model, [target]))
            (found if exists else absent).append(target)
        if found:
            report.matched[rule.name] = found
        if absent:
            report.missing[rule.name] = absent
        if rule.targets and not found:
            report.inapplicable.append(rule.name)
    return report
