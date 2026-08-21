"""A regulatory layer for genome-scale models: sensors → responses → constraint changes.

A stoichiometric model has no memory of *why* a reaction is available. Bounds are set once
and held, so a single model file can only ever represent one regulatory state. That is the
root of several problems already seen in this project:

* ``iSynCJ816_STAR`` ships its NDH-1 CO₂-uptake complexes bounded to zero — the
  **low-CO₂-inducible** systems, frozen in their repressed state. Forcing them open
  (as ``iSyn6803_STARplus`` does) freezes the opposite state. Neither is right across
  conditions; the correct object is a rule.
* The Calvin cycle runs in the dark, because nothing tells the model that
  ferredoxin–thioredoxin activation of PRK/GAPDH/FBPase/SBPase requires light.
* The biomass-specific absorption cross-section — the parameter the growth predictions
  are most sensitive to — is treated as a constant, when nitrogen starvation halves it by
  degrading the phycobilisomes.

**Design.** Each rule is a *sensor* (something readable from the environment or a flux), a
*response function* (how strongly the rule fires, 0–1), and a *target + effect*. Effects
reach bounds, enzyme costs, the protein budget and biomass composition, because
cyanobacterial regulation is mostly **resource reallocation** — and in an
enzyme-constrained model that is expressible rather than something to be caricatured as
an on/off switch.

Rules are evaluated to a fixed point: solve → read sensors → re-apply → solve, until the
state stops moving. Non-convergence is reported, never hidden.

**Two disciplines make this trustworthy rather than merely flexible:**

1. Every rule carries provenance and a confidence level, and every result records which
   rules fired and how hard. A prediction resting on an ``assumed`` threshold says so.
2. An empty rule set reproduces the unregulated model exactly. Regulation is opt-in.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import cobra

# --------------------------------------------------------------------------------------
# Confidence — how much weight a result resting on this rule deserves
# --------------------------------------------------------------------------------------
MEASURED = "measured"      # threshold and shape from primary literature
INFERRED = "inferred"      # mechanism established, quantities from related work
ASSUMED = "assumed"        # mechanism established, quantities chosen for plausibility

_CONFIDENCE_ORDER = {MEASURED: 0, INFERRED: 1, ASSUMED: 2}

#: Activation below which a rule is *not* considered to have fired.
#:
#: A saturating de-repression response never reaches exactly zero — NdhR at 17.6 mM
#: nitrate still evaluates to ~1e-4 — so testing ``activation > 0`` marks every such rule
#: as permanently active and hides the transitions that matter. 5 % is the point at which
#: a rule is doing something a reader would call regulation.
ACTIVE_THRESHOLD = 0.05


# --------------------------------------------------------------------------------------
# Response functions: sensor value → activation in [0, 1]
# --------------------------------------------------------------------------------------
def step(threshold: float, above: bool = True) -> Callable[[float], float]:
    """Hard switch. Use only where biology really is all-or-nothing (dark → CBB off)."""
    def f(x: float) -> float:
        return 1.0 if ((x >= threshold) == above) else 0.0
    return f


def hill(half: float, coefficient: float = 2.0, rising: bool = True
         ) -> Callable[[float], float]:
    """Saturating response — the usual shape for an induced or repressed system.

    ``half`` is the sensor value giving half-maximal response. ``rising=False`` gives
    de-repression as the sensor *falls*, which is how NdhR-controlled CCM genes behave.
    """
    def f(x: float) -> float:
        x = max(x, 0.0)
        if half <= 0:
            return 1.0 if rising else 0.0
        frac = (x ** coefficient) / (half ** coefficient + x ** coefficient)
        return frac if rising else 1.0 - frac
    return f


def ramp(low: float, high: float, rising: bool = True) -> Callable[[float], float]:
    """Linear between two points, clamped outside. Honest when the shape is unknown."""
    def f(x: float) -> float:
        if high == low:
            return 1.0 if x >= high else 0.0
        frac = (x - low) / (high - low)
        frac = min(1.0, max(0.0, frac))
        return frac if rising else 1.0 - frac
    return f


RESPONSES = {"step": step, "hill": hill, "ramp": ramp}


# --------------------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------------------
#: What a rule may change.
GATE = "gate"                # bounds → 0 when fired (or opened when not)
SCALE = "scale"              # capacity multiplied
ENZYME_COST = "enzyme_cost"  # prot_pool coefficient multiplied
BUDGET = "budget"            # total protein budget multiplied
BIOMASS = "biomass"          # biomass coefficient multiplied
PARAMETER = "parameter"      # a named simulation parameter (e.g. absorption cross-section)


@dataclass
class Rule:
    """One regulatory statement."""

    name: str
    regulator: str                       # the protein/system responsible, for the record
    sensor: str                          # key into the sensed state
    response: Callable[[float], float]
    effect: str                          # GATE | SCALE | ENZYME_COST | BUDGET | BIOMASS | PARAMETER
    targets: Sequence[str] = ()          # reaction ids, metabolite ids, or parameter names
    #: Value applied at full activation. For SCALE/ENZYME_COST/BUDGET/BIOMASS this is the
    #: multiplier at activation 1; the applied factor interpolates from 1.0 at activation 0.
    magnitude: float = 0.0
    description: str = ""
    provenance: str = ""
    confidence: str = ASSUMED
    enabled: bool = True
    #: The response as data (``{"kind": "hill", "half": 1.0, ...}``). Kept alongside the
    #: callable so a rule built in the editor round-trips to JSON without being rebuilt
    #: from a function object, which cannot be serialised.
    spec: dict = field(default_factory=dict)

    def factor(self, activation: float) -> float:
        """Multiplier at this activation, interpolating from 1.0 (off) to ``magnitude``."""
        return 1.0 + activation * (self.magnitude - 1.0)

    def describe_response(self) -> str:
        """One line saying when this rule fires, in the sensor's own units."""
        kind = self.spec.get("kind", "")
        if kind == "step":
            direction = "at or above" if self.spec.get("above", True) else "below"
            return f"fires fully {direction} {self.spec.get('threshold', 0):g}"
        if kind == "hill":
            rising = self.spec.get("rising", True)
            return (f"half-maximal at {self.spec.get('half', 0):g}, "
                    f"{'rising with' if rising else 'falling with'} the sensor "
                    f"(n={self.spec.get('coefficient', 2):g})")
        if kind == "ramp":
            low, high = self.spec.get("low", 0), self.spec.get("high", 0)
            if self.spec.get("rising", True):
                return f"ramps from 0 at {low:g} to full at {high:g}"
            return f"ramps from full at {low:g} down to 0 at {high:g}"
        return "response not described"


@dataclass
class RuleFiring:
    rule: str
    regulator: str
    sensor: str
    sensor_value: float
    activation: float
    effect: str
    applied: str
    confidence: str
    #: Did this rule actually change the model on this pass?
    #:
    #: Not the same as "activation is high". A gate is *closed* when activation is near
    #: zero — the CBB light rule shuts the Calvin cycle in the dark precisely because it
    #: is not activated. Judging by activation alone reported that rule as dormant while
    #: it was changing predicted growth by ~10%, which is the worst kind of silence.
    acted: bool = False

    def line(self) -> str:
        return (f"{self.rule} [{self.regulator}] — {self.sensor}={self.sensor_value:.4g} "
                f"→ {self.activation * 100:.0f}% · {self.applied} ({self.confidence})")


@dataclass
class RegulationResult:
    growth: float = 0.0
    status: str = ""
    firings: List[RuleFiring] = field(default_factory=list)
    parameters: Dict[str, float] = field(default_factory=dict)
    iterations: int = 0
    converged: bool = True
    notes: List[str] = field(default_factory=list)

    @staticmethod
    def is_active(firing: "RuleFiring") -> bool:
        """Whether a firing should be reported as having done something."""
        return bool(firing.acted) or firing.activation >= ACTIVE_THRESHOLD

    @property
    def weakest_confidence(self) -> str:
        active = [f.confidence for f in self.firings if RegulationResult.is_active(f)]
        if not active:
            return MEASURED
        return max(active, key=lambda c: _CONFIDENCE_ORDER.get(c, 9))

    def summary(self) -> str:
        if not self.firings:
            return "No regulatory rules fired — this is the unregulated model."
        fired = [f for f in self.firings if RegulationResult.is_active(f)]
        lines = [f"{len(fired)} of {len(self.firings)} rule(s) fired "
                 f"(converged in {self.iterations} iteration(s))."]
        lines += ["  • " + f.line() for f in fired]
        if self.weakest_confidence == ASSUMED:
            lines.append("  ⚠ At least one active rule rests on an ASSUMED threshold — "
                         "treat the quantitative result as indicative.")
        if not self.converged:
            lines.append("  ⚠ The regulatory state did not settle; the reported solution "
                         "is the last iterate, not a fixed point.")
        return "\n".join(lines)


class RuleSet:
    """A collection of rules, applied together.

    ``organism`` is carried deliberately. Regulation is organism-specific in a way
    stoichiometry is not, and a set that does not say what it describes invites being
    applied to the wrong model.
    """

    def __init__(self, rules: Optional[Sequence[Rule]] = None, name: str = "",
                 organism: str = "", description: str = "", source: str = ""):
        self.rules: List[Rule] = list(rules or [])
        self.name = name
        self.organism = organism
        self.description = description
        self.source = source

    def __len__(self) -> int:
        return len(self.rules)

    def enabled(self) -> List[Rule]:
        return [r for r in self.rules if r.enabled]

    def by_regulator(self) -> Dict[str, List[Rule]]:
        out: Dict[str, List[Rule]] = {}
        for rule in self.rules:
            out.setdefault(rule.regulator, []).append(rule)
        return out


# --------------------------------------------------------------------------------------
# Sensors
# --------------------------------------------------------------------------------------
#: The sensors a rule can read that mean something beyond one exchange reaction.
#:
#: A cell does not sense "EX_no3_e"; it senses nitrogen availability, whichever source
#: supplies it. Keeping these as named sensors is what lets one NtcA rule stay correct
#: when the user switches the culture from nitrate to ammonium to urea.
STANDARD_SENSORS: Dict[str, str] = {
    "light_uE": "Light (µmol photons m⁻² s⁻¹)",
    "inorganic_c_mM": "Inorganic carbon (mM)",
    "nitrogen_mM": "Nitrogen, any source (mM)",
    "iron_uM": "Iron (µM)",
    "phosphate_mM": "Phosphate (mM)",
    "sulfur_mM": "Sulfur (mM)",
    "growth": "Growth rate (h⁻¹, from the previous solve)",
    "c_to_n_flux": "Carbon:nitrogen uptake flux ratio",
    "photon_flux": "Photon uptake flux (mmol gDW⁻¹ h⁻¹)",
}

#: Metabolite-id fragments that identify which standard sensor an exchange feeds.
#: Order matters: the first match wins, so specific ids precede generic ones.
_SENSOR_BY_METABOLITE = [
    (("photon", "hnu", "light"), "light_uE"),
    (("co2", "hco3", "co3"), "inorganic_c_mM"),
    (("no3", "no2", "nh4", "nh3", "urea", "gln__l", "glu__l", "n2"), "nitrogen_mM"),
    (("fe2", "fe3", "fe"), "iron_uM"),
    (("pi", "po4", "phosphate"), "phosphate_mM"),
    (("so4", "so3", "h2s", "cys__l"), "sulfur_mM"),
]


def exchange_metabolite_id(exchange_id: str) -> str:
    """The metabolite an exchange id refers to, e.g. ``EX_no3_e`` → ``no3_e``."""
    text = exchange_id
    for prefix in ("EX_", "EX", "DM_", "SK_"):
        if text.startswith(prefix):
            text = text[len(prefix):]
            break
    return text


def sensor_for_exchange(exchange_id: str) -> str:
    """Which sensor an exchange should drive.

    Falls back to a per-metabolite key (``no3_e_mM``) when the compound is not one of the
    standard nutrient classes. That fallback is deliberately *not* silent-nothing: a rule
    can still be written against it, it simply has no organism-level meaning attached.
    """
    met = exchange_metabolite_id(exchange_id).lower()
    stem = met[:-2] if met.endswith("_e") else met
    for fragments, sensor in _SENSOR_BY_METABOLITE:
        if stem in fragments:
            return sensor
    for fragments, sensor in _SENSOR_BY_METABOLITE:
        if any(stem.startswith(f) for f in fragments):
            return sensor
    return f"{met}_mM"


def sensor_label(sensor: str) -> str:
    return STANDARD_SENSORS.get(sensor, f"{sensor} (custom)")


# --------------------------------------------------------------------------------------
# Sensing
# --------------------------------------------------------------------------------------
def sense(model: cobra.Model, environment: Dict[str, float],
          solution=None) -> Dict[str, float]:
    """Build the sensed state a rule set reads.

    Environment values are what the caller set (light, available N, available Ci).
    Flux-derived values need a solution and therefore an iteration; they are omitted on
    the first pass, which is why the fixed point exists.
    """
    state = dict(environment)
    if solution is not None and getattr(solution, "status", "") == "optimal":
        def flux(rid: str) -> float:
            return float(solution.fluxes.get(rid, 0.0) or 0.0)
        state["growth"] = float(solution.objective_value or 0.0)
        # A carbon:nitrogen flux ratio is the closest stoichiometric proxy for the
        # 2-oxoglutarate signal that NtcA/PII actually read.
        carbon = abs(flux("EX_co2_e")) + abs(flux("EX_hco3_e")) + abs(flux("EX_glc__D_e"))
        nitrogen = abs(flux("EX_no3_e")) + abs(flux("EX_nh4_e"))
        state["c_to_n_flux"] = carbon / nitrogen if nitrogen > 1e-9 else float("inf")
        state["photon_flux"] = abs(flux("EX_photon_e"))
    return state


# --------------------------------------------------------------------------------------
# Application
# --------------------------------------------------------------------------------------
def _targets_in(model: cobra.Model, ids: Sequence[str]) -> List[cobra.Reaction]:
    found = []
    for rid in ids:
        if model.reactions.has_id(rid):
            found.append(model.reactions.get_by_id(rid))
        else:
            # Tolerate AutoPACMEN decoration (FOO → FOO_1 / FOO_TG_forward).
            import re
            base = re.sub(r"_TG_(forward|reverse)$", "", rid)
            base = re.sub(r"_\d+$", "", base)
            found.extend(r for r in model.reactions
                         if re.sub(r"_\d+$", "",
                                   re.sub(r"_TG_(forward|reverse)$", "", r.id)) == base)
    return found


def apply_rules(model: cobra.Model, ruleset: RuleSet, state: Dict[str, float],
                parameters: Dict[str, float]) -> List[RuleFiring]:
    """Apply every enabled rule to ``model`` in place. Returns what fired."""
    from .flux_context import FLUX_TOL  # noqa: F401 - shared tolerance convention

    firings: List[RuleFiring] = []
    prot_pool = (model.metabolites.get_by_id("prot_pool")
                 if model.metabolites.has_id("prot_pool") else None)

    for rule in ruleset.enabled():
        if rule.sensor not in state:
            continue                     # sensor not available on this pass
        value = state[rule.sensor]
        if value == float("inf"):
            activation = 1.0
        else:
            activation = max(0.0, min(1.0, rule.response(value)))
        factor = rule.factor(activation)
        applied = ""
        acted = False

        if rule.effect == GATE:
            reactions = _targets_in(model, rule.targets)
            # activation 1 = fully ON; the rule states what the regulator permits.
            for rxn in reactions:
                if activation < 0.5:
                    rxn.lower_bound = min(0.0, max(rxn.lower_bound, 0.0))
                    rxn.upper_bound = 0.0
                    if rxn.lower_bound < 0:
                        rxn.lower_bound = 0.0
            closed = activation < 0.5 and bool(reactions)
            applied = (f"{len(reactions)} reaction(s) "
                       + ("open" if activation >= 0.5 else "CLOSED"))
            acted = closed

        elif rule.effect == SCALE:
            reactions = _targets_in(model, rule.targets)
            for rxn in reactions:
                if rxn.upper_bound > 0:
                    rxn.upper_bound *= factor
                if rxn.lower_bound < 0:
                    rxn.lower_bound *= factor
            applied = f"{len(reactions)} reaction(s) × {factor:.3g}"
            acted = bool(reactions) and abs(factor - 1.0) > 0.01

        elif rule.effect == ENZYME_COST:
            if prot_pool is None:
                continue
            reactions = _targets_in(model, rule.targets)
            changed = 0
            for rxn in reactions:
                if prot_pool in rxn.metabolites:
                    current = rxn.metabolites[prot_pool]
                    # add_metabolites is additive: pass the delta needed to reach
                    # current × factor.
                    rxn.add_metabolites({prot_pool: current * (factor - 1.0)})
                    changed += 1
            applied = f"{changed} enzyme cost(s) × {factor:.3g}"
            acted = changed > 0 and abs(factor - 1.0) > 0.01

        elif rule.effect == BUDGET:
            if model.reactions.has_id("ER_pool_TG_"):
                pool = model.reactions.get_by_id("ER_pool_TG_")
                pool.upper_bound *= factor
                applied = f"protein budget × {factor:.3g} → {pool.upper_bound:.4g}"
                acted = abs(factor - 1.0) > 0.01

        elif rule.effect == BIOMASS:
            changed = 0
            for rxn in model.reactions:
                if "BIOMASS" not in rxn.id.upper():
                    continue
                for mid in rule.targets:
                    if model.metabolites.has_id(mid):
                        met = model.metabolites.get_by_id(mid)
                        if met in rxn.metabolites:
                            current = rxn.metabolites[met]
                            rxn.add_metabolites({met: current * (factor - 1.0)})
                            changed += 1
            applied = f"{changed} biomass coefficient(s) × {factor:.3g}"
            acted = changed > 0 and abs(factor - 1.0) > 0.01

        elif rule.effect == PARAMETER:
            for pname in rule.targets:
                parameters[pname] = parameters.get(pname, 1.0) * factor
            applied = f"{', '.join(rule.targets)} × {factor:.3g}"
            acted = abs(factor - 1.0) > 0.01

        firings.append(RuleFiring(
            rule=rule.name, regulator=rule.regulator, sensor=rule.sensor,
            sensor_value=value, activation=activation, effect=rule.effect,
            applied=applied, confidence=rule.confidence, acted=acted))
    return firings


def simulate(model: cobra.Model, ruleset: RuleSet, environment: Dict[str, float], *,
             max_iterations: int = 8, tolerance: float = 1e-5) -> RegulationResult:
    """Solve with regulation applied, iterating to a fixed point.

    With an empty rule set this is a plain FBA solve and returns exactly what the
    unregulated model would.
    """
    result = RegulationResult()
    previous = None
    solution = None

    for iteration in range(1, max_iterations + 1):
        with model as scratch:
            parameters: Dict[str, float] = {}
            state = sense(scratch, environment, solution)
            firings = apply_rules(scratch, ruleset, state, parameters)
            solution = scratch.optimize()
            growth = float(solution.objective_value or 0.0)
            result.status = str(solution.status)
            result.growth = growth
            result.firings = firings
            result.parameters = parameters
            result.iterations = iteration

        if previous is not None and abs(growth - previous) < tolerance:
            result.converged = True
            break
        previous = growth
    else:
        result.converged = len(ruleset) == 0
        if not result.converged:
            result.notes.append(
                f"did not settle within {max_iterations} iterations "
                f"(last change {abs((previous or 0) - result.growth):.3g})")
    return result


# --------------------------------------------------------------------------------------
# Serialisation — rules live in an editable file, not in code
# --------------------------------------------------------------------------------------
def to_dict(ruleset: RuleSet) -> dict:
    out = {
        "name": ruleset.name,
        "organism": getattr(ruleset, "organism", ""),
        "description": getattr(ruleset, "description", ""),
        "source": getattr(ruleset, "source", ""),
        "rules": [],
    }
    for rule in ruleset.rules:
        out["rules"].append({
            "name": rule.name, "regulator": rule.regulator, "sensor": rule.sensor,
            "effect": rule.effect, "targets": list(rule.targets),
            "magnitude": rule.magnitude, "description": rule.description,
            "provenance": rule.provenance, "confidence": rule.confidence,
            "enabled": rule.enabled,
            "response": dict(rule.spec or {}),
        })
    return out


def build_response(spec: dict) -> Callable[[float], float]:
    """Turn a response spec into the callable a rule evaluates."""
    kind = spec.get("kind", "step")
    return RESPONSES[kind](**{k: v for k, v in spec.items() if k != "kind"})


def from_dict(payload: dict) -> RuleSet:
    rules = []
    for item in payload.get("rules", []):
        spec = item.get("response") or {}
        rule = Rule(
            name=item["name"], regulator=item.get("regulator", ""),
            sensor=item["sensor"], response=build_response(spec),
            effect=item["effect"],
            targets=item.get("targets", []), magnitude=item.get("magnitude", 0.0),
            description=item.get("description", ""),
            provenance=item.get("provenance", ""),
            confidence=item.get("confidence", ASSUMED),
            enabled=item.get("enabled", True),
            spec=dict(spec))
        rules.append(rule)
    return RuleSet(rules, name=payload.get("name", ""),
                   organism=payload.get("organism", ""),
                   description=payload.get("description", ""),
                   source=payload.get("source", ""))


def load(path: str) -> RuleSet:
    with open(path, encoding="utf-8") as fh:
        return from_dict(json.load(fh))


def save(ruleset: RuleSet, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(to_dict(ruleset), fh, indent=2)
