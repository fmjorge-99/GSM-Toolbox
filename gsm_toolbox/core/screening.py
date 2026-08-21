"""Screen a model across conditions, and follow a culture through time.

Steady-state FBA answers "what can the cell do *here*". Most regulation, though, is about
what the cell does as conditions *change* — nitrogen running out over a batch, inorganic
carbon falling as a dense culture draws it down, light shifting between day and night.
This module provides the two shapes of experiment that make those questions askable:

**Scan** — sweep one or two environmental variables over a grid, solve at each point, and
report growth, target fluxes and which regulatory rules fired. Answers "how does the
network rewire across a nitrogen gradient".

**Timecourse** — a dynamic FBA: integrate biomass and substrate concentrations forward,
re-solving at each step with the regulatory state re-evaluated from the *current*
concentrations. Answers "what happens over a batch as nitrogen is consumed", which is the
question a steady-state solve cannot reach.

Both take the regulatory rule set, so a scan shows regulation switching as conditions
cross a threshold rather than a smooth stoichiometric response.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import cobra

from . import regulation as reg


# --------------------------------------------------------------------------------------
# Condition scanning
# --------------------------------------------------------------------------------------
@dataclass
class ScanPoint:
    environment: Dict[str, float]
    growth: float = 0.0
    status: str = ""
    fluxes: Dict[str, float] = field(default_factory=dict)
    active_rules: List[str] = field(default_factory=list)
    weakest_confidence: str = reg.MEASURED
    protein_saturation: Optional[float] = None


@dataclass
class ScanResult:
    variables: List[str]
    points: List[ScanPoint] = field(default_factory=list)
    targets: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def table(self) -> List[dict]:
        rows = []
        for p in self.points:
            row = {v: p.environment.get(v) for v in self.variables}
            row["growth"] = round(p.growth, 6)
            for t in self.targets:
                row[t] = round(p.fluxes.get(t, 0.0), 6)
            row["rules_fired"] = len(p.active_rules)
            row["confidence"] = p.weakest_confidence
            if p.protein_saturation is not None:
                row["protein_saturation"] = round(p.protein_saturation, 4)
            rows.append(row)
        return rows

    def transitions(self, key: str = "growth", threshold: float = 0.15) -> List[dict]:
        """Points where the response changes sharply — where regulation switches.

        A scan is most informative at its breakpoints; this finds them rather than
        leaving the user to read a table of numbers.
        """
        out = []
        rows = self.table()
        for a, b in zip(rows, rows[1:]):
            before, after = a.get(key) or 0.0, b.get(key) or 0.0
            if before <= 1e-9:
                continue
            change = abs(after - before) / before
            if change >= threshold:
                out.append({"from": {v: a[v] for v in self.variables},
                            "to": {v: b[v] for v in self.variables},
                            key: (before, after),
                            "relative_change": round(change, 3)})
        return out


def _apply_environment(model: cobra.Model, environment: Dict[str, float],
                       mapping: Dict[str, Tuple[str, float]]) -> None:
    """Translate environment values into exchange bounds via ``mapping``.

    ``mapping`` maps an environment key to (exchange id, conversion factor from the
    environment unit to mmol gDW⁻¹ h⁻¹).
    """
    for key, (exchange, factor) in mapping.items():
        if key not in environment or not model.reactions.has_id(exchange):
            continue
        model.reactions.get_by_id(exchange).lower_bound = -abs(environment[key] * factor)


def scan(model: cobra.Model, ruleset: reg.RuleSet, *,
         variables: Dict[str, Sequence[float]],
         base_environment: Dict[str, float],
         mapping: Dict[str, Tuple[str, float]],
         targets: Sequence[str] = (),
         objective: Optional[str] = None) -> ScanResult:
    """Sweep one or two environment variables and record the response at each point."""
    names = list(variables)
    grids = [list(variables[n]) for n in names]
    result = ScanResult(variables=names, targets=list(targets))

    def combos(index: int, current: Dict[str, float]):
        if index == len(names):
            yield dict(current)
            return
        for value in grids[index]:
            current[names[index]] = value
            yield from combos(index + 1, current)

    for environment_values in combos(0, {}):
        environment = dict(base_environment)
        environment.update(environment_values)
        point = ScanPoint(environment=dict(environment_values))
        with model as scratch:
            if objective and scratch.reactions.has_id(objective):
                scratch.objective = objective
            _apply_environment(scratch, environment, mapping)
            outcome = reg.simulate(scratch, ruleset, environment)
            point.growth = outcome.growth
            point.status = outcome.status
            point.active_rules = [f.rule for f in outcome.firings
                                  if reg.RegulationResult.is_active(f)]
            point.weakest_confidence = outcome.weakest_confidence
            # Re-solve once with the settled regulatory state to read target fluxes.
            if outcome.status == "optimal" and targets:
                state = reg.sense(scratch, environment, None)
                reg.apply_rules(scratch, ruleset, state, {})
                solution = scratch.optimize()
                for t in targets:
                    if scratch.reactions.has_id(t):
                        point.fluxes[t] = float(solution.fluxes.get(t, 0.0) or 0.0)
                if scratch.reactions.has_id("ER_pool_TG_"):
                    pool = scratch.reactions.get_by_id("ER_pool_TG_")
                    used = abs(float(solution.fluxes.get("ER_pool_TG_", 0.0) or 0.0))
                    if pool.upper_bound:
                        point.protein_saturation = used / pool.upper_bound
        result.points.append(point)
    return result


# --------------------------------------------------------------------------------------
# Dynamic simulation (batch / nutrient depletion)
# --------------------------------------------------------------------------------------
@dataclass
class Substrate:
    """A medium component that is consumed, and how it maps into the model.

    ``sensor`` is what regulation reads and is kept separate from ``label`` on purpose.
    A culture switched from nitrate to ammonium is still nitrogen-limited, so both
    exchanges drive the ``nitrogen_mM`` sensor and one NtcA rule stays correct across the
    switch. Without that indirection every rule would have to name an exchange, and a rule
    set would silently stop applying the moment the user changed the nitrogen source.
    """

    name: str                 # environment key the rules sense (the sensor)
    exchange: str             # exchange reaction id
    initial_mM: float
    #: mmol gDW⁻¹ h⁻¹ of uptake per mM present — the affinity of the uptake system.
    #: A saturating (Monod) form is used, so this is the maximum specific uptake.
    max_uptake: float
    half_saturation_mM: float = 0.05
    label: str = ""           # human name for the table column
    #: Supplied in excess: uptake keeps the medium's own bound and the concentration is
    #: held constant. For a component that is buffered or replenished — bicarbonate in a
    #: buffered medium, a nutrient deliberately kept non-limiting. Without this mode the
    #: only way to include such a component would be to let it deplete under a guessed
    #: affinity, which invents a limitation the experiment does not have.
    buffered: bool = False

    def column(self) -> str:
        return self.label or self.name

    def uptake_limit(self, concentration: float) -> float:
        if concentration <= 0:
            return 0.0
        return self.max_uptake * concentration / (self.half_saturation_mM + concentration)


#: Above this, an exchange bound is a modelling placeholder for "unlimited", not a
#: measured capacity — cobra writes 1000 for anything left open.
_UNBOUNDED = 100.0

#: Specific uptake used when the model states no real capacity. Of the order of a
#: measured maximal uptake rate (glucose in *E. coli* is ~10 mmol gDW⁻¹ h⁻¹).
DEFAULT_MAX_UPTAKE = 10.0


def default_uptake(model: cobra.Model, exchange: str) -> float:
    """A usable specific uptake for this exchange, in mmol gDW⁻¹ h⁻¹.

    Steady-state models leave uptake bounds at 1000 to mean "not limiting". Carried into
    a dynamic run that is not a capacity but an instruction to consume the entire medium
    within one time step, which produces a spike and then a dead culture. Treating the
    placeholder as a placeholder — and substituting a physiological rate the user can see
    and change — is what keeps the time course interpretable.
    """
    bound = 0.0
    if model.reactions.has_id(exchange):
        bound = abs(model.reactions.get_by_id(exchange).lower_bound)
    if not bound or bound >= _UNBOUNDED:
        return DEFAULT_MAX_UPTAKE
    return bound


def substrate_from_exchange(model: cobra.Model, exchange: str, initial_mM: float, *,
                            sensor: str = "", max_uptake: float = 0.0,
                            buffered: bool = False) -> Substrate:
    """Build a :class:`Substrate` for any exchange the model has.

    ``max_uptake`` defaults to the capacity the medium grants, unless that is cobra's
    unlimited placeholder — see :func:`default_uptake`. The regulatory sensor is derived
    from the compound unless overridden.
    """
    from . import regulation as reg

    if not max_uptake:
        max_uptake = default_uptake(model, exchange)
    return Substrate(
        name=sensor or reg.sensor_for_exchange(exchange),
        exchange=exchange,
        initial_mM=float(initial_mM),
        max_uptake=float(max_uptake),
        label=exchange_label(model, exchange),
        buffered=buffered)


def exchange_label(model: cobra.Model, exchange: str) -> str:
    """A readable name for an exchange, or its id when the name adds nothing.

    Several published models store the id itself as the name, sometimes with the
    parentheses SBML-escaped — iJN678 calls `EX_hco3_e` "EX hco3 LPAREN e RPAREN".
    Tidying that string still leaves an id pretending to be a name, so the id is used
    directly instead; it is what the reader recognises anyway.
    """
    from . import regulation as reg
    from .physiology import clean_label

    if not model.reactions.has_id(exchange):
        return reg.exchange_metabolite_id(exchange)
    name = clean_label(model.reactions.get_by_id(exchange).name or "").strip()
    if not name:
        return exchange
    # "EX hco3(e)" against "EX_hco3_e": the same characters once separators are dropped.
    squashed = "".join(ch for ch in name.lower() if ch.isalnum())
    if squashed == "".join(ch for ch in exchange.lower() if ch.isalnum()):
        return exchange
    return name


@dataclass
class Phase:
    """What the culture experiences at one moment — sensors and bounds together.

    The two travel in one object on purpose. `environment` is what a regulatory rule
    perceives; `uptake_limits` is what the solver is physically allowed. Supplying one
    without the other silently compares two different experiments: a darkness the cell
    senses but the network does not, or a photon bound the rules never learn about.
    """
    #: Sensor name → value, merged over the run's base environment.
    environment: Dict[str, float] = field(default_factory=dict)
    #: Exchange id → maximum uptake rate (a magnitude; the sign is applied here).
    uptake_limits: Dict[str, float] = field(default_factory=dict)


def diel_schedule(*, day: Phase, night: Phase, day_h: float = 12.0,
                  night_h: float = 12.0, start_in_day: bool = True
                  ) -> Callable[[float], Phase]:
    """A repeating light/dark cycle, the condition a cyanobacterium actually grows in.

    Returned as a callable of time so `timecourse` stays agnostic about *why* conditions
    change; any other driver (a feed, a temperature ramp) has the same shape.
    """
    period = float(day_h) + float(night_h)
    if period <= 0:
        raise ValueError("a diel cycle needs a positive day or night length")

    def phase_at(time_h: float) -> Phase:
        into = math.fmod(float(time_h), period)
        in_day = into < day_h
        return (day if in_day else night) if start_in_day else (night if in_day else day)

    return phase_at


@dataclass
class TimePoint:
    time_h: float
    biomass_gDW_L: float
    concentrations: Dict[str, float]      # column label → mM
    growth_rate: float
    fluxes: Dict[str, float] = field(default_factory=dict)
    active_rules: List[str] = field(default_factory=list)
    #: Secretion rate of the product exchange, mmol gDW-1 h-1 (0 when none is requested).
    product_rate: float = 0.0
    #: Product accumulated in the medium up to this point, mM.
    product_mM: float = 0.0


@dataclass
class TimecourseResult:
    points: List[TimePoint] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    #: Column names by role, so a plot can offer sensible axes without re-deriving them.
    concentration_columns: List[str] = field(default_factory=list)
    flux_columns: List[str] = field(default_factory=list)
    #: Exchange whose secretion was followed, when the run asked for one.
    product: Optional[str] = None

    def table(self) -> List[dict]:
        rows = []
        for p in self.points:
            row = {"time_h": round(p.time_h, 3),
                   "biomass_gDW_L": round(p.biomass_gDW_L, 5),
                   "growth_rate_per_h": round(p.growth_rate, 5)}
            row.update({f"{k} (mM)": round(v, 5) for k, v in p.concentrations.items()})
            row.update({k: round(v, 5) for k, v in p.fluxes.items()})
            if self.product:
                row[f"{self.product} rate"] = round(p.product_rate, 5)
                row[f"{self.product} (mM)"] = round(p.product_mM, 5)
            row["rules_fired"] = ";".join(p.active_rules)
            rows.append(row)
        return rows

    def numeric_columns(self) -> List[str]:
        """Columns a plot can use for an axis — everything except the rule listing."""
        rows = self.table()
        if not rows:
            return []
        return [c for c in rows[0] if c != "rules_fired"]

    def phase_changes(self) -> List[dict]:
        """Times at which the set of active rules changes — regulatory transitions."""
        out = []
        for a, b in zip(self.points, self.points[1:]):
            before, after = set(a.active_rules), set(b.active_rules)
            if before != after:
                out.append({
                    "time_h": round(b.time_h, 2),
                    "gained": sorted(after - before),
                    "lost": sorted(before - after),
                    "growth_before": round(a.growth_rate, 5),
                    "growth_after": round(b.growth_rate, 5)})
        return out


def timecourse(model: cobra.Model, ruleset: reg.RuleSet, *,
               substrates: Sequence[Substrate],
               base_environment: Dict[str, float],
               initial_biomass: float = 0.05,
               duration_h: float = 120.0,
               step_h: float = 1.0,
               targets: Sequence[str] = (),
               objective: Optional[str] = None,
               schedule: Optional[Callable[[float], Phase]] = None,
               product: Optional[str] = None,
               product_growth_fraction: float = 0.9) -> TimecourseResult:
    """Dynamic FBA: follow a batch culture as it consumes its substrates.

    At each step the substrate concentrations set the uptake bounds *and* the regulatory
    sensors, so a rule fires when the culture crosses its threshold — which is the point
    of doing this dynamically rather than at steady state.

    `schedule` drives conditions that change with the clock rather than with consumption
    — a light/dark cycle above all (see `diel_schedule`). It returns a `Phase`, which
    carries the sensor values and the uptake bounds together so the two cannot drift
    apart.

    `product` names an exchange whose secretion is followed and accumulated. Growth alone
    never selects a product an engineered strain was built to make, so the step is solved
    twice: maximise growth, hold it at `product_growth_fraction` of that maximum, then
    maximise the product. The fraction is the engineering assumption — how much growth the
    strain is allowed to give up — and it is recorded in the notes rather than buried.

    Integration is explicit Euler on biomass and concentrations. That is adequate for the
    qualitative question ("when does the culture switch state") and is *not* adequate for
    precise kinetics; the step size is reported so the reader can judge.
    """
    result = TimecourseResult()
    result.notes.append(
        f"Explicit Euler, step {step_h} h — suitable for phase behaviour, not for "
        "precise kinetics.")
    if product:
        result.notes.append(
            f"Product {product} maximised at {product_growth_fraction:.0%} of the "
            "growth rate achievable at each step.")
    result.concentration_columns = [f"{s.column()} (mM)" for s in substrates]
    result.flux_columns = list(targets)
    result.product = product if product and model.reactions.has_id(product) else None

    biomass = initial_biomass
    # Keyed by exchange, not by sensor: two nitrogen sources are two pools that deplete
    # independently even though they drive the same regulatory sensor.
    concentrations = {s.exchange: float(s.initial_mM) for s in substrates}
    time = 0.0
    product_mM = 0.0

    def sensed(current: Dict[str, float]) -> Dict[str, float]:
        """Sensor values from the current pools — several sources summing into one.

        Nitrate plus ammonium is not two half-starvations; the cell reads total available
        nitrogen. Summing per sensor is what makes a mixed-source culture behave sensibly.
        """
        totals: Dict[str, float] = {}
        for s in substrates:
            totals[s.name] = totals.get(s.name, 0.0) + current[s.exchange]
        return totals

    while time <= duration_h + 1e-9:
        phase = schedule(time) if schedule else None
        environment = dict(base_environment)
        environment.update(sensed(concentrations))
        if phase is not None:
            # The clock overrides the pools: a component held at a fixed concentration
            # still goes dark at night.
            environment.update(phase.environment)

        with model as scratch:
            if objective and scratch.reactions.has_id(objective):
                scratch.objective = objective
            # Uptake capacity follows a Monod form of the remaining concentration.
            # A buffered component keeps whatever the medium already allowed.
            for s in substrates:
                if s.buffered or not scratch.reactions.has_id(s.exchange):
                    continue
                limit = s.uptake_limit(concentrations[s.exchange])
                scratch.reactions.get_by_id(s.exchange).lower_bound = -limit
            if phase is not None:
                for rid, limit in phase.uptake_limits.items():
                    if scratch.reactions.has_id(rid):
                        scratch.reactions.get_by_id(rid).lower_bound = -abs(limit)

            outcome = reg.simulate(scratch, ruleset, environment)
            growth = outcome.growth if outcome.status == "optimal" else 0.0

            fluxes: Dict[str, float] = {}
            uptakes: Dict[str, float] = {}
            product_rate = 0.0
            if outcome.status == "optimal":
                state = reg.sense(scratch, environment, None)
                reg.apply_rules(scratch, ruleset, state, {})
                # A second, product-directed solve on the already-regulated network.
                # Regulation is applied before the product objective, not after, so a
                # gate the cell has closed cannot be reopened to serve the product.
                if product and scratch.reactions.has_id(product):
                    growth_reaction = _objective_reaction(scratch)
                    if growth_reaction is not None:
                        growth_reaction.lower_bound = max(
                            growth_reaction.lower_bound,
                            product_growth_fraction * growth)
                    scratch.objective = product
                solution = scratch.optimize()
                if solution.status != "optimal":
                    growth = 0.0
                else:
                    if product and scratch.reactions.has_id(product):
                        product_rate = float(solution.fluxes.get(product, 0.0) or 0.0)
                    for t in targets:
                        if scratch.reactions.has_id(t):
                            fluxes[t] = float(solution.fluxes.get(t, 0.0) or 0.0)
                    for s in substrates:
                        uptakes[s.exchange] = abs(float(
                            solution.fluxes.get(s.exchange, 0.0) or 0.0))

        result.points.append(TimePoint(
            time_h=time, biomass_gDW_L=biomass,
            concentrations={s.column(): concentrations[s.exchange] for s in substrates},
            growth_rate=growth,
            fluxes=fluxes,
            active_rules=[f.rule for f in outcome.firings
                          if reg.RegulationResult.is_active(f)],
            product_rate=product_rate,
            product_mM=product_mM))

        # Integrate. Substrate consumed = specific uptake × biomass × dt.
        # A buffered pool is held: it is being replenished as fast as it is used, which
        # is what "supplied in excess" means. Draining it while leaving uptake unbounded
        # would report a concentration the simulation is not actually respecting.
        for s in substrates:
            if s.buffered:
                continue
            consumed = uptakes.get(s.exchange, 0.0) * biomass * step_h
            concentrations[s.exchange] = max(
                0.0, concentrations[s.exchange] - consumed)
        product_mM += product_rate * biomass * step_h
        biomass *= math.exp(growth * step_h)
        time += step_h

    return result


def _objective_reaction(model: cobra.Model) -> Optional["cobra.Reaction"]:
    """The single reaction the objective maximises, if it is a single reaction.

    A lexicographic product solve has to hold growth at a floor, which means finding the
    reaction that carries it. A composite objective has no one reaction to pin, so the
    caller is told nothing was found rather than being given the wrong one.
    """
    carrying = [r for r in model.reactions if r.objective_coefficient]
    return carrying[0] if len(carrying) == 1 else None
