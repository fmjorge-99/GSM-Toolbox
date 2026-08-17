"""Thermodynamic feasibility of a designed pathway: per-reaction ΔrG′ and the
Max-min Driving Force (MDF).

A route that is stoichiometrically perfect can still be thermodynamically hopeless: if
every metabolite concentration that makes step A go forward makes step B go backward,
no steady state carries flux in the intended direction. The MDF (Noor et al., PLoS
Comput. Biol. 2014) is the standard single number for this — it optimises the metabolite
concentrations to make the *least* favourable reaction as favourable as possible:

    ΔrG′_j(c) = ΔrG′°_j + RT · Σ_i S_ij · ln c_i
    MDF       = max_c  min_j  ( −ΔrG′_j(c) )      subject to  c_min ≤ c_i ≤ c_max

A positive MDF (kJ/mol) means concentrations exist at which every step is downhill —
the pathway is feasible; MDF ≤ 0 means at least one step is stuck no matter the
concentrations. The optimisation is a linear program in {ln c_i} and the objective B.

The LP here is self-contained and dependency-light (SciPy only). The hard part —
obtaining ΔrG′° for each reaction — is delegated to eQuilibrator when available; that
package and its compound cache are optional and downloaded on demand (see
``equilibrator_available`` / ``reaction_dg0``).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

R_KJ = 8.314e-3        # gas constant, kJ / (mol·K)
T_DEFAULT = 298.15     # K
DEFAULT_CONC = (1e-6, 1e-2)   # physiological metabolite range, molar (1 µM … 10 mM)


class ThermoDataError(RuntimeError):
    """The thermodynamics engine or its data could not be loaded.

    Distinct from "this route has no ΔrG′°": that is a normal result about the chemistry,
    whereas this means the analysis could not run at all and the message should tell the
    user how to repair their installation.
    """
REFERENCE_CONC = 1e-3         # 1 mM — the conventional reference for reporting ΔrG′


def _dg_at_reference(dg0: float, stoich: dict, rt: float,
                     fixed: Optional[dict] = None) -> float:
    """ΔrG′ with every free metabolite at :data:`REFERENCE_CONC` (1 mM).

    This is the meaningful way to report a *single* reaction: unlike an MDF it does not
    optimise concentrations, so two reactions can be compared directly.
    """
    import math
    fixed = fixed or {}
    total = dg0
    for met, coeff in stoich.items():
        conc = fixed.get(met, REFERENCE_CONC)
        if conc and conc > 0:
            total += rt * coeff * math.log(conc)
    return total


@dataclass
class MDFResult:
    """Outcome of a Max-min Driving Force analysis over a pathway."""

    mdf: float = float("nan")                 # kJ/mol; > 0 ⇒ feasible
    feasible: bool = False
    bottleneck: str = ""                      # reaction with the least driving force
    dg_prime: Dict[str, float] = field(default_factory=dict)   # ΔrG′ at the MDF optimum
    dg0_prime: Dict[str, float] = field(default_factory=dict)  # standard ΔrG′°
    ln_conc: Dict[str, float] = field(default_factory=dict)    # optimised ln concentrations
    missing: List[str] = field(default_factory=list)  # reactions with no ΔrG′° available
    note: str = ""
    single_reaction: bool = False   # magnitude not comparable to a multi-step MDF (L3)
    # For a single reaction: ΔrG′ with every free metabolite at 1 mM. This, not `mdf`, is
    # the number to show and compare in that case (VI.7).
    dg_prime_reference: Optional[float] = None
    #: reaction id → why it could not be scored, for a warning that names the gap
    missing_reasons: Dict[str, str] = field(default_factory=dict)
    #: ``(scored, total)`` reactions — how much of the route the number actually covers
    coverage: Tuple[int, int] = (0, 0)

    @property
    def coverage_fraction(self) -> float:
        scored, total = self.coverage
        return (scored / total) if total else 0.0

    def data_warning(self) -> str:
        """Plain statement of what the database could not supply, or "" when complete.

        Deliberately names the reactions and the compounds behind them: "insufficient
        data" that does not say *which* data leaves the user with nowhere to go.
        """
        scored, total = self.coverage
        if not total or scored == total:
            return ""
        lines = [f"<b>{total - scored} of {total} reaction(s) in this route have no "
                 "thermodynamic data</b> in the loaded databases."]
        if scored == 0:
            lines.append("Nothing could be scored, so no driving force can be computed "
                         "for this route at all.")
        else:
            lines.append(f"The driving force below covers only the {scored} reaction(s) "
                         "that could be scored — treat it as a partial answer.")
        detail = [f"&nbsp;&nbsp;• <b>{rid}</b> — {why}"
                  for rid, why in list(self.missing_reasons.items())[:8]]
        if detail:
            lines.append("<br>".join(detail))
        if len(self.missing_reasons) > 8:
            lines.append(f"…and {len(self.missing_reasons) - 8} more.")
        return "<br><br>".join(lines)


def max_min_driving_force(
        stoich: Dict[str, Dict[str, float]],
        dg0: Dict[str, float],
        *,
        conc_bounds: Optional[Dict[str, Tuple[float, float]]] = None,
        default_conc: Tuple[float, float] = DEFAULT_CONC,
        fixed: Optional[Dict[str, float]] = None,
        temperature: float = T_DEFAULT) -> MDFResult:
    """Solve the MDF linear program.

    ``stoich`` maps reaction id → {metabolite id: coefficient} (products positive,
    substrates negative), already oriented in the intended *forward* direction. ``dg0``
    maps reaction id → ΔrG′° (kJ/mol). ``conc_bounds`` optionally overrides the molar
    concentration range per metabolite; ``fixed`` pins metabolites (e.g. water, or a
    held cofactor ratio) to an exact molar concentration. Reactions absent from ``dg0``
    are ignored for the objective and reported in ``missing``.
    """
    import numpy as np
    from scipy.optimize import linprog

    rxns = [r for r in stoich if r in dg0]
    missing = [r for r in stoich if r not in dg0]
    res = MDFResult(dg0_prime={r: dg0[r] for r in rxns}, missing=missing)
    if not rxns:
        res.note = "No ΔrG′° available for any reaction in this route."
        return res

    mets = sorted({m for r in rxns for m in stoich[r]})
    idx = {m: i for i, m in enumerate(mets)}
    rt = R_KJ * temperature
    fixed = fixed or {}
    conc_bounds = conc_bounds or {}

    # Variables: x = [ln c_0 … ln c_{n-1}, B].  Maximise B  ⇔  minimise −B.
    n = len(mets)
    c = np.zeros(n + 1)
    c[-1] = -1.0

    # For each reaction:  −ΔrG′_j ≥ B
    #   −(ΔG0_j + RT Σ S_ij ln c_i) ≥ B
    #   RT Σ S_ij ln c_i + B ≤ −ΔG0_j
    A_ub, b_ub = [], []
    for r in rxns:
        row = np.zeros(n + 1)
        for m, s in stoich[r].items():
            row[idx[m]] = rt * s
        row[-1] = 1.0
        A_ub.append(row)
        b_ub.append(-dg0[r])

    bounds = []
    for m in mets:
        if m in fixed and fixed[m] > 0:
            lc = float(np.log(fixed[m]))
            bounds.append((lc, lc))
        else:
            lo, hi = conc_bounds.get(m, default_conc)
            bounds.append((float(np.log(lo)), float(np.log(hi))))
    bounds.append((None, None))    # B is free

    sol = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds,
                  method="highs")
    if not sol.success:
        res.note = f"MDF optimisation did not converge ({sol.message})."
        return res

    ln_c = {m: float(sol.x[idx[m]]) for m in mets}
    res.ln_conc = ln_c
    res.mdf = float(sol.x[-1])
    res.feasible = res.mdf > 1e-9
    # A one-reaction "pathway" has no shared intermediate to trade off, so the LP simply
    # pins every substrate at its ceiling and the product at its floor. The resulting
    # number can run to hundreds of kJ/mol and is NOT comparable to a multi-step MDF.
    # Flag it so the caller can present it honestly rather than as a headline value.
    res.single_reaction = len(rxns) == 1
    if res.single_reaction:
        # An MDF over one reaction is not a pathway driving force at all: with no shared
        # intermediate to trade off, the optimiser simply pins every substrate at its
        # ceiling and the product at its floor, producing values in the hundreds. Report
        # ΔrG′ at REFERENCE concentrations instead — a number that means something and is
        # comparable between single reactions (VI.7).
        rid = rxns[0]
        res.dg_prime_reference = _dg_at_reference(dg0[rid], stoich[rid], rt, fixed)
        res.note = ("Single reaction, not a pathway: an MDF needs shared intermediates to "
                    "trade off, so the optimiser here simply pins substrates high and the "
                    "product low and the magnitude is meaningless. Use the ΔrG′ at "
                    f"reference concentrations ({REFERENCE_CONC * 1e3:g} mM) instead.")
    # ΔrG′ of each reaction at the optimum, and the tightest (bottleneck) one.
    worst, worst_df = "", float("inf")
    for r in rxns:
        dg = dg0[r] + rt * sum(s * ln_c[m] for m, s in stoich[r].items())
        res.dg_prime[r] = dg
        if -dg < worst_df:
            worst_df, worst = -dg, r
    res.bottleneck = worst
    return res


# --------------------------------------------------------------------------------------
# ΔrG′° sourcing via eQuilibrator (optional, downloaded on demand)
# --------------------------------------------------------------------------------------

# Metabolite annotation namespaces eQuilibrator can resolve, best identifier first.
_ID_NAMESPACES = [
    ("kegg.compound", "kegg"),
    ("bigg.metabolite", "bigg.metabolite"),
    ("metanetx.chemical", "metanetx.chemical"),
    ("seed.compound", "seed.compound"),
    ("chebi", "CHEBI"),
    ("inchi_key", "inchikey"),
]


CACHE_BYTES = 1_340_000_000       # ~1.34 GB — the eQuilibrator compound cache
CACHE_DOI = "10.5281/zenodo.4128543"


def mdf_suite_enabled() -> bool:
    """True only if the user has opted into the thermodynamics suite in Preferences.

    The MDF feature depends on a ~1.34 GB external dataset and on assumptions (pH, ionic
    strength) a user must consciously choose, so it is hidden by default. Every MDF entry
    point in the GUI is gated on this.
    """
    from . import preferences
    return preferences.mdf_enabled()


def mdf_ready() -> bool:
    """Enabled AND actually usable (package importable and compound cache present)."""
    return mdf_suite_enabled() and equilibrator_available() and cache_present()


def cache_path() -> Optional[str]:
    """Where eQuilibrator keeps (or would keep) its compound cache."""
    try:
        import pooch
        return os.path.join(str(pooch.os_cache("equilibrator")), "compounds.sqlite")
    except Exception:  # noqa: BLE001
        return None


def cache_present() -> bool:
    p = cache_path()
    if p and os.path.exists(p) and os.path.getsize(p) > 1_000_000:
        return True
    return _bundled_cache_path() is not None


def download_equilibrator_cache(progress=None) -> str:
    """Download the eQuilibrator compound cache (~1.34 GB) on explicit user consent.

    ``progress`` is an optional ``callable(fraction, message)``. Returns the cache path.
    Raises RuntimeError with a user-facing message on failure. This is only ever called
    from the Preferences opt-in flow — never implicitly.
    """
    if not equilibrator_available():
        raise RuntimeError(
            "The eQuilibrator package is not installed in this build, so the "
            "thermodynamics suite cannot be enabled.")
    ensure_equilibrator_cache()          # a bundled copy, if any, satisfies this
    if cache_present():
        return cache_path() or ""
    if progress:
        progress(0.0, "Downloading the eQuilibrator compound database (~1.34 GB)…")
    try:
        # Constructing ComponentContribution triggers the pooch-managed download of
        # compounds.sqlite from Zenodo, with resume/checksum handling.
        _component_contribution()
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "Could not download the eQuilibrator compound database. Check your internet "
            f"connection and try again.\n\n{exc}") from exc
    if progress:
        progress(1.0, "Thermodynamics data ready.")
    return cache_path() or ""


def equilibrator_available() -> bool:
    """True if the eQuilibrator API can be imported (the compound cache may still need
    to download on first real use)."""
    import importlib.util
    return importlib.util.find_spec("equilibrator_api") is not None


def parameters_present() -> bool:
    """True when the component-contribution parameter file is already downloaded.

    Separate from :func:`cache_present`: the compound database ships with the app, but
    the ~69 MB parameter file does not and is fetched on first use. Knowing which is
    missing lets the UI say "this is a one-time download" instead of appearing to hang.
    """
    try:
        import pooch
        path = os.path.join(str(pooch.os_cache("equilibrator")), "cc_params.npz")
        return os.path.exists(path) and os.path.getsize(path) > 1_000_000
    except Exception:  # noqa: BLE001
        return False


def _bundled_cache_path() -> Optional[str]:
    """Path to a ``compounds.sqlite`` shipped inside the app, if present.

    PyInstaller unpacks bundled data under ``sys._MEIPASS``; in a source checkout it may
    sit next to the package. Returns None if no bundled cache is found."""
    import sys
    candidates = []
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        candidates.append(os.path.join(meipass, "equilibrator_data", "compounds.sqlite"))
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates.append(os.path.join(here, "resources", "equilibrator_data",
                                   "compounds.sqlite"))
    for c in candidates:
        if os.path.exists(c) and os.path.getsize(c) > 1_000_000:
            return c
    return None


def ensure_equilibrator_cache() -> None:
    """Put the bundled compound cache where eQuilibrator expects it, avoiding the ~1.3 GB
    runtime download when the app ships it. No-op if eQuilibrator already has the cache,
    or if nothing is bundled (it then downloads on demand as usual)."""
    bundled = _bundled_cache_path()
    if bundled is None:
        return
    try:
        import pooch
        # pooch.os_cache("equilibrator") already resolves to the …/equilibrator/Cache dir
        # where equilibrator-cache expects compounds.sqlite directly.
        cache_dir = str(pooch.os_cache("equilibrator"))
        dest = os.path.join(cache_dir, "compounds.sqlite")
        if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
            return
        os.makedirs(cache_dir, exist_ok=True)
        import shutil
        shutil.copy2(bundled, dest)
    except Exception:  # noqa: BLE001 — fall back to eQuilibrator's own download
        pass


def _accessions(met) -> List[str]:
    """Every eQuilibrator-resolvable accession for a metabolite, best first.

    Returns a *list*, and reads annotation keys in any spelling. Both matter:

    * the old version looked only for MIRIAM-style keys (``kegg.compound``,
      ``inchi_key``), so on a BiGG-derived database — which writes ``"KEGG Compound"``
      and ``"InChI Key"``, holding identifiers.org URLs — it resolved **nothing**. Every
      compound then fell through to a per-compound online structure lookup, which is why
      the analysis was slow and scored so little;
    * it also returned only the first match, so one failed lookup lost a compound that a
      second identifier would have resolved.
    """
    from . import identifiers

    found = identifiers.normalised(met)
    out: List[str] = []
    # InChIKey first: it identifies a structure exactly, so it cannot resolve to the
    # wrong compound the way a database id occasionally can.
    for key in found.get("inchikey", []):
        out.append(f"inchikey:{key}")
    for value in found.get("kegg", []):
        out.append(f"kegg:{value}")
    for value in found.get("metanetx", []):
        out.append(f"metanetx.chemical:{value}")
    for value in found.get("chebi", []):
        out.append(value if value.upper().startswith("CHEBI") else f"CHEBI:{value}")
    for value in found.get("bigg", []):
        out.append(f"bigg.metabolite:{value}")
    for value in found.get("seed", []):
        out.append(f"seed.compound:{value}")
    for value in found.get("hmdb", []):
        out.append(f"hmdb:{value}")
    if not out:
        # A bare BiGG id is still worth a try when nothing is annotated at all.
        base = (met.id or "").rsplit("_", 1)[0]
        if base:
            out.append(f"bigg.metabolite:{base}")
    return out


def _accession(met) -> Optional[str]:
    """First resolvable accession, kept for callers that want a single value."""
    found = _accessions(met)
    return found[0] if found else None


_CC = None      # ComponentContribution is expensive (~65 s); build once per process.


def _component_contribution():
    """The cached predictor. Raises ``ThermoDataError`` with an actionable message.

    Constructing this downloads and unpickles eQuilibrator's parameter file, and a
    failure here aborts the whole analysis — so the error must say what actually went
    wrong rather than surfacing a library-internal message like "No module named
    'numpy._core'", which tells the user nothing they can act on.
    """
    global _CC
    if _CC is None:
        from . import numpy_compat
        # Its parameter file is a numpy 2-era pickle; under the numpy 1 this app pins
        # that needs the compatibility aliases in place first.
        numpy_compat.install()
        try:
            from equilibrator_api import ComponentContribution
            _CC = ComponentContribution()
        except ModuleNotFoundError as exc:
            if "numpy" in str(exc):
                raise ThermoDataError(
                    "This build of the app cannot read eQuilibrator's thermodynamic "
                    "parameter file: a NumPy compatibility module is missing from the "
                    f"package ({exc}). Reinstalling the latest version fixes it."
                ) from exc
            raise ThermoDataError(
                f"A component of the thermodynamics engine is missing: {exc}") from exc
        except Exception as exc:  # noqa: BLE001
            raise ThermoDataError(
                "The thermodynamics data could not be loaded — it may be missing or "
                f"only partly downloaded. Details: {_short(exc, 160)}") from exc
    return _CC


def reaction_dg0(model, reaction_ids: List[str], *,
                 temperature: float = T_DEFAULT,
                 ph: float = 7.5, ionic_strength_M: float = 0.25,
                 pmg: float = 3.0,
                 offline: bool = False) -> Tuple[Dict[str, float], List[str]]:
    """Standard transformed reaction energies ΔrG′° (kJ/mol) via eQuilibrator.

    Returns ``(dg0_by_reaction, unresolved_reaction_ids)``. Raises ``RuntimeError`` with
    a user-facing message if eQuilibrator is not installed. The first call constructs the
    component-contribution predictor, which downloads the compound cache (~hundreds of MB)
    once — hence this must run off the UI thread.
    """
    if not equilibrator_available():
        raise RuntimeError(
            "Thermodynamic analysis needs the optional eQuilibrator add-on, which is "
            "not installed.")
    ensure_equilibrator_cache()     # use the bundled cache if the app ships one
    from equilibrator_api import Q_

    cc = _component_contribution()   # cached: constructed once per process (~65 s)
    cc.p_h = Q_(ph)
    cc.ionic_strength = Q_(f"{ionic_strength_M} M")
    cc.p_mg = Q_(pmg)
    cc.temperature = Q_(f"{temperature} K")

    # Resolve every metabolite once (compound cache lookups are the expensive part).
    # `diagnosis` records WHY anything failed, so the caller can tell the user which
    # compounds lack data instead of reporting a bare "no thermodynamics available".
    compounds: Dict[str, object] = {}
    diagnosis: Dict[str, str] = {}
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            continue
        for met in model.reactions.get_by_id(rid).metabolites:
            if met.id in compounds:
                continue
            cpd, why = _resolve_compound(cc, met, offline=offline)
            compounds[met.id] = cpd
            if cpd is None:
                diagnosis[met.id] = why

    dg0: Dict[str, float] = {}
    unresolved: List[str] = []
    reasons: Dict[str, str] = {}
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            unresolved.append(rid)
            reasons[rid] = "not present in the model being analysed"
            continue
        rxn = model.reactions.get_by_id(rid)
        parts = []
        blockers = []
        for met, coeff in rxn.metabolites.items():
            cpd = compounds.get(met.id)
            if cpd is None:
                blockers.append(met.name or met.id)
                continue
            parts.append((coeff, cpd))
        if blockers:
            # Name the compounds, not just the reaction: that is what the user has to act
            # on (add an identifier, fetch a structure, or accept the gap).
            unresolved.append(rid)
            reasons[rid] = ("no thermodynamic data for " + ", ".join(blockers[:4])
                            + (" …" if len(blockers) > 4 else ""))
            continue
        try:
            eq_rxn = _phased_reaction(cc, parts)
            dg = cc.standard_dg_prime(eq_rxn)
            dg0[rid] = float(dg.value.m_as("kJ/mol"))
        except Exception as exc:  # noqa: BLE001
            unresolved.append(rid)
            reasons[rid] = f"eQuilibrator could not score it ({_short(exc)})"
    _LAST_REASONS.clear()
    _LAST_REASONS.update(reasons)
    return dg0, unresolved


#: Why each reaction of the most recent call could not be scored; read by the GUI.
_LAST_REASONS: Dict[str, str] = {}


def last_failure_reasons() -> Dict[str, str]:
    """Per-reaction explanations from the most recent :func:`reaction_dg0` call."""
    return dict(_LAST_REASONS)


def _short(exc: Exception, limit: int = 80) -> str:
    text = str(exc).strip().replace("\n", " ") or type(exc).__name__
    return text[:limit] + ("…" if len(text) > limit else "")


def _resolve_compound(cc, met, *, offline: bool = False):
    """``(compound, reason)`` — try every identifier, then the structure.

    Returns the first accession that resolves. Trying them all matters: a compound whose
    KEGG id is missing from the cache is often present under its InChIKey or ChEBI id,
    and the previous single-accession attempt threw those away.
    """
    tried = 0
    for acc in _accessions(met):
        tried += 1
        try:
            cpd = cc.get_compound(acc)
        except Exception:  # noqa: BLE001 - unknown accession, malformed id, cache miss
            cpd = None
        if cpd is not None:
            return cpd, ""
    # No identifier resolved. eQuilibrator can also decompose a compound from its
    # STRUCTURE, which is the only route for rule-generated intermediates (VI.5).
    cpd = _compound_from_structure(cc, met, offline=offline)
    if cpd is not None:
        return cpd, ""
    if tried:
        return None, (f"none of its {tried} identifier(s) are in the eQuilibrator "
                      "cache, and its structure could not be decomposed")
    return None, "it has no usable identifier and no structure"


def _compound_from_structure(cc, met, *, offline: bool = False):
    """Resolve a compound from its 2-D structure when no accession works (VI.5).

    eQuilibrator's component-contribution method decomposes a molecule into groups, so it
    can estimate ΔfG′° for a compound it has never seen — provided it is given a
    structure. Rule-generated intermediates always have SMILES, and database metabolites
    can have one fetched, so this is what restores thermodynamics for novel chemistry.

    ``offline`` skips the web lookup and uses only structures already annotated or cached.
    This is now a *fallback* rather than the main path: with identifiers read correctly
    almost everything resolves from the local cache, so the slow per-compound web request
    is reserved for genuinely novel chemistry.
    """
    try:
        from .chemistry import metabolite_smiles
        smi = metabolite_smiles(met, online=not offline)
        if not smi:
            return None
    except Exception:  # noqa: BLE001
        return None
    for getter, arg in (("get_compound_by_inchi", _inchi_from_smiles(smi)),
                        ("get_compound", f"smiles:{smi}")):
        if not arg:
            continue
        fn = getattr(cc, getter, None)
        if fn is None:
            continue
        try:
            cpd = fn(arg)
            if cpd is not None:
                return cpd
        except Exception:  # noqa: BLE001
            continue
    return None


def _inchi_from_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchi(m) if m is not None else ""
    except Exception:  # noqa: BLE001
        return ""


def _phased_reaction(cc, parts):
    """Build an eQuilibrator Reaction from (coefficient, compound) pairs."""
    from equilibrator_api import Reaction
    stoich = {}
    for coeff, cpd in parts:
        stoich[cpd] = stoich.get(cpd, 0) + coeff
    return Reaction(stoich)


def analyse_pathway_mdf(model, reaction_ids: List[str], *,
                        conc_bounds: Optional[Dict[str, Tuple[float, float]]] = None,
                        default_conc: Tuple[float, float] = DEFAULT_CONC,
                        **dg0_kwargs) -> MDFResult:
    """End-to-end MDF for a route in ``model``: source ΔrG′°, then solve the LP.

    Water and protons are held at unit activity (they are part of the transformed ΔrG′°,
    not free concentration variables). ``default_conc`` is the allowed concentration range
    (M) for metabolites without an explicit entry in ``conc_bounds``; widening it makes the
    MDF more optimistic, so it is a user-visible assumption rather than a constant.
    Raises ``RuntimeError`` if eQuilibrator is absent.
    """
    dg0, unresolved = reaction_dg0(model, reaction_ids, **dg0_kwargs)
    stoich: Dict[str, Dict[str, float]] = {}
    fixed: Dict[str, float] = {}
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            continue
        rxn = model.reactions.get_by_id(rid)
        row = {}
        for met, coeff in rxn.metabolites.items():
            base = (met.name or met.id).lower()
            if base.startswith(("h2o", "water")) or met.id.lower().startswith("h2o") \
                    or base in ("h+", "proton") or met.id.lower().startswith("h_"):
                fixed[met.id] = 1.0            # activity folded into ΔrG′°
            row[met.id] = coeff
        stoich[rid] = row
    res = max_min_driving_force(stoich, dg0, conc_bounds=conc_bounds,
                                default_conc=default_conc, fixed=fixed,
                                temperature=dg0_kwargs.get("temperature", T_DEFAULT))
    res.missing = unresolved
    res.missing_reasons = last_failure_reasons()
    res.coverage = (len(dg0), len(dg0) + len(unresolved))
    return res
