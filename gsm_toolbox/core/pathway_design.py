"""De novo / heterologous pathway design (OptStrain-style).

Rule-based retrosynthesis suites (PathPred, BNICE.ch, RetroPath) rely on large
external reaction-rule and compound-structure databases. The tractable equivalent
inside a genome-scale model is to draw candidate reactions from a *universal*
reaction database (e.g. the BiGG universal model or MetaNetX) and find the minimal
set of heterologous reactions that lets the host produce a target metabolite —
the first step of the OptStrain workflow.

This module predicts such pathways (via COBRApy gap-filling against the universal
database) and can apply the chosen reactions to the working model so the rest of
the toolbox (FBA, strain design, …) can be used on the engineered host.
"""

from __future__ import annotations

import os
import re
import urllib.request
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import cobra
import pandas as pd
from cobra.flux_analysis import gapfill

# BiGG Models provides a COBRA-compatible universal model (all reactions/metabolites)
# and per-model JSON files. See http://bigg.ucsd.edu/data_access
BIGG_UNIVERSAL_URL = "http://bigg.ucsd.edu/static/namespace/universal_model.json"
BIGG_MODEL_URL = "http://bigg.ucsd.edu/static/models/{model_id}.json"


class PathwayDesignError(Exception):
    """Raised when a heterologous pathway cannot be predicted or applied."""


class DatabaseDownloadError(Exception):
    """Raised when an online reaction database cannot be downloaded."""


def _cache_dir() -> str:
    base = os.path.join(os.path.expanduser("~"), ".gsm_toolbox", "databases")
    os.makedirs(base, exist_ok=True)
    return base


def _download(url: str, dest: str) -> str:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "GSM-ToolBox"})
        with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as fh:
            fh.write(resp.read())
    except Exception as exc:  # noqa: BLE001
        raise DatabaseDownloadError(f"Download failed for {url}:\n{exc}") from exc
    return dest


def download_bigg_universal(force: bool = False) -> cobra.Model:
    """Download (and cache) the BiGG universal reaction database, returning a model."""
    from .io_models import load_model

    dest = os.path.join(_cache_dir(), "bigg_universal_model.json")
    if force or not os.path.exists(dest):
        _download(BIGG_UNIVERSAL_URL, dest)
    return load_model(dest)


def bundled_universal_available() -> bool:
    """Whether a curated offline universal ships with this install (Issue 4)."""
    from .. import resources
    return resources.has_offline_universal()


def load_bundled_universal() -> cobra.Model:
    """Load the curated offline universal that ships with the app (no network).

    This is the always-available fallback when online databases can't be reached.
    """
    from .. import resources
    from .io_models import load_model
    path = resources.offline_universal_path()
    if not os.path.exists(path):
        raise DatabaseDownloadError(
            "No offline universal is bundled with this build. Load a reaction "
            "database from a file instead, or connect to the internet to fetch one.")
    return load_model(path)


def fetch_bigg_model(model_id: str) -> cobra.Model:
    """Download a specific BiGG model by id (e.g. 'iML1515') to use as a reaction source."""
    from .io_models import load_model

    model_id = model_id.strip()
    if not model_id:
        raise DatabaseDownloadError("Provide a BiGG model id (e.g. iML1515).")
    dest = os.path.join(_cache_dir(), f"bigg_{model_id}.json")
    if not os.path.exists(dest):
        _download(BIGG_MODEL_URL.format(model_id=model_id), dest)
    return load_model(dest)


@dataclass
class PathwayResult:
    target: str
    reactions: pd.DataFrame = field(default_factory=pd.DataFrame)  # id, name, equation, source
    n_steps: int = 0
    production_flux: float = float("nan")
    note: str = ""
    reaction_ids: List[str] = field(default_factory=list)
    # Correctness disclosures (Issue 2): human-readable warnings about any
    # mass/charge-unbalanced or unverifiable steps in this route.
    warnings: List[str] = field(default_factory=list)
    balanced: bool = True
    # Database ids of the step(s) that introduce a dead-end metabolite and so stop this
    # route carrying flux. Used to retry the search without the offending step while
    # keeping the rest of the route (forbidding the whole route would also throw away
    # the essential terminal enzyme it shares with the working alternative).
    blocked_by: List[str] = field(default_factory=list)
    # When no route was found: the compound(s) nothing in the database can make. These
    # are exactly what an online fetch has to supply, so the GUI can offer to go and
    # get the missing chemistry rather than leaving the user stuck.
    missing_compounds: List[str] = field(default_factory=list)
    # True for rule-derived (RetroRules) routes. Their `production_flux` is measured on a
    # synthetic demand for a novel compound whose only link to metabolism is the proposed
    # route itself, so it saturates against a generic host bound rather than expressing a
    # validated pathway capacity — unrelated products returned the identical figure. The
    # GUI must therefore present it qualitatively ("carries flux" / "no flux"), never as a
    # number comparable with a database route (L2).
    flux_is_indicative: bool = False
    # Carbon yield (mol C in product / mol C consumed), when it could be computed — a
    # normalised figure that IS comparable across targets, unlike raw flux (L10).
    carbon_yield: float = float("nan")
    # Why `carbon_yield` is NaN, when it is. NaN alone reads as "zero yield" or "bug";
    # a missing product formula is neither, and the user can act on it.
    carbon_yield_note: str = ""
    # Steps whose mass/charge balance could not be verified (a missing formula, not an
    # imbalance). Reported separately so the headline warning stays meaningful (VI.14a).
    unverified_steps: List[str] = field(default_factory=list)
    # step id -> the participants whose missing formula blocked its balance check. Naming
    # them is what makes "unverified" actionable: it is the shopping list for a fetch.
    unverified_reasons: Dict[str, List[str]] = field(default_factory=dict)
    # Skeletal-isomer check (VI.1A): human-readable descriptions of steps whose carbon
    # backbone changes in a way no single reaction can achieve. Empty means either clean
    # or not yet checked — `isomer_checked` distinguishes the two.
    isomer_warnings: List[str] = field(default_factory=list)
    isomer_checked: bool = False
    isomer_coverage: str = ""            # e.g. "2/3" reactions actually checkable
    # Database reaction ids collapsed into this route because they are the same chemistry
    # with the same cofactors (VI.12). Kept so the user can still see the alternatives.
    equivalent_routes: List[str] = field(default_factory=list)

    def duplicate(self) -> "PathwayResult":
        """An independent copy of this route, for editing as a separate variant.

        The common case is one core pathway that needs several different terminal steps:
        copy it, add a different final reaction to each copy, and compare them. The copy
        must share nothing mutable with the original — including the reactions table and
        the cached analyses hung off it — or editing one variant would silently change
        the other.
        """
        import copy as _copy

        clone = _copy.deepcopy(self)
        # Analyses (diagnosis, branching, MDF) describe the route as it was; a variant
        # that is about to be edited must be re-analysed rather than inherit them.
        for attr in ("_diagnosis", "_branching", "_mdf"):
            if hasattr(clone, attr):
                delattr(clone, attr)
        return clone

    def carbon_yield_text(self) -> str:
        """The carbon yield as a phrase, or why it could not be computed.

        Never render ``carbon_yield`` raw: NaN in a yield column reads as a failed
        design, when the usual cause is a database metabolite with no formula.
        """
        cy = self.carbon_yield
        if cy == cy:                                   # not NaN
            return f"{cy * 100:.1f}% of consumed C"
        return f"not computable — {self.carbon_yield_note}" if self.carbon_yield_note \
            else "not computable"

    def missing_formula_labels(self) -> List[str]:
        """Every participant, across all steps, whose missing formula blocks a check.

        ``'name (id)'`` labels, in first-seen order — exactly what an on-demand formula
        fetch has to supply for this route.
        """
        out: List[str] = []
        for labels in self.unverified_reasons.values():
            for label in labels:
                if label not in out:
                    out.append(label)
        return out


_BALANCE_TOL = 1e-6


def reaction_balance_verdict(rxn: cobra.Reaction):
    """The full :class:`~gsm_toolbox.core.balancing.BalanceVerdict` for ``rxn``.

    Use this rather than :func:`reaction_balance` whenever the answer is shown to a
    user: it names the participant whose missing formula made the check impossible,
    which is the difference between "this route loses 42 carbons" and "we have no
    formula for its product".
    """
    # Formula-aware (back-fills formulas from annotation) so far fewer reactions
    # are reported as "unverifiable" than a raw formula-attribute check would give.
    from . import balancing
    return balancing.assess_reaction(rxn)


def reaction_balance(rxn: cobra.Reaction):
    """Assess a reaction's mass/charge balance.

    Returns ``(balanced, residual, checkable)`` where ``residual`` maps each
    unbalanced element (or ``'charge'``) to its non-zero net coefficient, and
    ``checkable`` is True only when every metabolite carries a formula (so the
    verdict is trustworthy). Reactions we cannot check are neither trusted nor
    dropped — they are flagged as "unverified" for the user.
    """
    return reaction_balance_verdict(rxn).as_tuple()


def _is_grossly_unbalanced(rxn: cobra.Reaction) -> bool:
    """True when a reaction is verifiably (formulas present) mass/charge-unbalanced.

    Used to keep physically impossible shortcuts (e.g. a C3 → C20 one-step jump)
    out of the pathway search entirely."""
    return reaction_balance_verdict(rxn).verifiably_imbalanced


def _format_residual(residual: dict) -> str:
    return ", ".join(f"{k}{v:+.3g}" for k, v in sorted(residual.items()))


def _ensure_metabolite(work: cobra.Model, universal: cobra.Model, met_id: str) -> cobra.Metabolite:
    """Return the target metabolite in ``work``, copying it from ``universal`` if absent."""
    if work.metabolites.has_id(met_id):
        return work.metabolites.get_by_id(met_id)
    if not universal.metabolites.has_id(met_id):
        raise PathwayDesignError(f"Metabolite '{met_id}' is in neither the model nor the database.")
    src = universal.metabolites.get_by_id(met_id)
    new = cobra.Metabolite(src.id, name=src.name, formula=src.formula,
                           charge=src.charge, compartment=src.compartment)
    work.add_metabolites([new])
    return work.metabolites.get_by_id(met_id)


def _base_key(met_id: str) -> str:
    """A compartment-independent key for a metabolite: its id without the trailing
    ``_<compartment>`` suffix. Used so a reaction is treated as the same chemistry
    regardless of which compartment the database placed it in."""
    return met_id.rsplit("_", 1)[0] if "_" in met_id else met_id


def _is_currency_id(met_id: str) -> bool:
    from .network_graph import CURRENCY_BASES
    return _base_key(met_id).lower() in CURRENCY_BASES


def _is_currency_key(key: str) -> bool:
    from .network_graph import CURRENCY_BASES
    return key.lower() in CURRENCY_BASES


# Cofactor / currency metabolites recognised by NAME, so detection works across
# namespaces (BiGG ids, MetaNetX MNXM ids with no formula, KEGG C-ids…). Matched
# after normalising away charges/isotopes/stereo. NOTE: acyl-CoA thioesters
# (malonyl-CoA, coumaroyl-CoA…) are carbon skeletons, NOT currency — only bare
# "CoA"/"coenzyme A" is — so the acyl-CoA guard below is essential.
_CURRENCY_NAMES = frozenset({
    "h2o", "water", "h", "h+", "proton", "hydron", "oh",
    "co2", "carbon dioxide", "o2", "dioxygen", "oxygen", "hydrogen peroxide",
    "coa", "coenzyme a",
    "atp", "adp", "amp", "adenosine monophosphate", "adenosine diphosphate",
    "adenosine triphosphate", "gtp", "gdp", "gmp", "ctp", "cdp", "cmp",
    "utp", "udp", "ump", "itp", "idp", "imp",
    "nad", "nadh", "nad+", "nadp", "nadph", "nadp+", "fad", "fadh2", "fmn", "fmnh2",
    "pi", "phosphate", "orthophosphate", "hydrogenphosphate",
    "ppi", "diphosphate", "pyrophosphate",
    "nh3", "nh4", "nh4+", "ammonia", "ammonium",
    "so4", "sulfate", "so3", "sulfite", "hco3", "bicarbonate", "co3", "carbonate",
    "h2s", "hydrogen sulfide", "hs", "no3", "nitrate", "no2", "nitrite",
    # --- generic redox carriers / electron shuttles -------------------------
    # These are regenerated by the host and must be treated as freely available,
    # otherwise whole classes of chemistry never "fire" during the pathway search.
    # In particular P450 / flavin-monooxygenase steps — exactly how terpenes are
    # oxidised (valencene -> nootkatol -> nootkatone) — are written in ModelSEED as
    # "O2 + Reduced_flavin + substrate <=> Flavin + H2O + H + product". Without these
    # entries "Reduced flavin" looked like a carbon skeleton that had to be
    # synthesised, so the reaction could never fire and the target was unreachable.
    "flavin", "reduced flavin", "oxidized flavin", "flavin reduced", "flavin oxidized",
    "fadh", "fmnh", "riboflavin",
    "acceptor", "reduced acceptor", "oxidized acceptor", "electron acceptor",
    "donor", "reduced donor", "oxidized donor", "electron donor",
    "ferredoxin", "reduced ferredoxin", "oxidized ferredoxin",
    "flavodoxin", "reduced flavodoxin", "oxidized flavodoxin",
    "thioredoxin", "reduced thioredoxin", "oxidized thioredoxin",
    "glutaredoxin", "glutathione", "reduced glutathione", "oxidized glutathione",
    "cytochrome c", "ferricytochrome c", "ferrocytochrome c",
    "cytochrome c oxidized", "cytochrome c reduced",
    "plastoquinone", "plastoquinol", "ubiquinone", "ubiquinol",
    "menaquinone", "menaquinol", "quinone", "quinol",
    "plastocyanin", "reduced plastocyanin", "oxidized plastocyanin",
    "photon", "light", "hnu",
    "s adenosyl l methionine", "s adenosyl l homocysteine", "sam", "sah",
    # --- generic CLASS placeholders --------------------------------------------
    # Some databases write a reaction against a whole class of compound rather than a
    # named one. ModelSEED's "NAD-P-OR-NOP"/"NADH-P-OR-NOP" (a stand-in for NAD(P)+ /
    # NAD(P)H) appears in 363 reactions EACH: unrecognised, every one of them was
    # unusable, because the search believed a compound named "NAD-P-OR-NOP" had to be
    # synthesised from scratch first.
    "nad p or nop", "nadh p or nop", "nadp or nop", "nadph or nop",
    "amino acids", "l amino acids", "amino acids 20", "d amino acids",
    "a protein", "protein", "holo acp", "holo acyl carrier protein",
})

# Substring patterns for generic redox carriers whose names vary a lot between
# databases (e.g. "Reduced acceptor", "Oxidized ferredoxin [iron-sulfur] cluster").
# Applied only to names that are NOT acyl-conjugates (guarded in _is_currency_met).
_CURRENCY_PATTERNS = (
    "ferredoxin", "flavodoxin", "thioredoxin", "glutaredoxin", "plastocyanin",
    "cytochrome", "plastoquino", "ubiquino", "menaquino",
    # Generic class placeholders (see the note in _CURRENCY_NAMES). Matched as
    # substrings because the exact wording varies a lot between databases.
    "or-nop", "or_nop",
    # Protein-bound residues ("Protein-Histidines", "a [protein]-L-lysine",
    # "Protein N(pi)-phospho-L-histidine"): the protein is part of the enzyme system
    # and is always present, so these are carriers rather than compounds to synthesise.
    "[protein]", "protein-histidin", "protein histidin", "protein-l-", "protein l-",
    "protein-n-", "protein n(", "protein lysine", "protein-lysine",
    "amino-acid residue", "amino acid residue",
)


def _currency_name_candidates(name: str) -> set:
    """Several normalised forms of a metabolite name for currency matching, so we can
    tolerate charges/isotopes/locants: e.g. ``coenzyme A``, ``ammonia-(13)N``,
    ``adenosine 5'-monophosphate``, ``((13)C)carbon dioxide``, ``H(+)``."""
    s = (name or "").strip().lower()
    prev = None
    while prev != s:                                 # peel nested parentheticals
        prev = s
        s = re.sub(r"\([^()]*\)", "", s)
    s = re.sub(r"[’'`]", "", s)
    base = re.sub(r"\s+", " ", re.sub(r"[-–,]+", " ", s)).strip()
    cands = {s.strip(), base}
    cands.add(re.sub(r"\s[a-z]$", "", base).strip())          # drop trailing isotope letter
    cands.add(re.sub(r"\s+", " ", re.sub(r"\d+", " ", base)).strip())  # drop locant digits
    return {c for c in cands if c}


def _is_currency_met(met) -> bool:
    """Whether a metabolite is a currency/cofactor — by id base OR by name, so it
    works even for formula-less MetaNetX/KEGG metabolites (fixes the scheme drawing
    treating H2O/CoA/CO2/ATP as pathway intermediates) and for generic redox carriers
    like ModelSEED's "Reduced flavin" (without which P450/monooxygenase steps never
    fire during the pathway search)."""
    if _is_currency_id(getattr(met, "id", "") or ""):
        return True
    raw = (getattr(met, "name", "") or "").strip().lower()
    if not raw:
        return False
    # Acyl-CoA / acyl-ACP thioesters are real carbon skeletons; only the bare carrier
    # is currency. Guard these BEFORE the generic patterns below.
    if "coa" in raw and raw not in ("coa", "coenzyme a"):
        return False
    if "acp" in raw or "acyl-carrier" in raw or "acyl carrier" in raw:
        return False
    cands = _currency_name_candidates(raw)
    if cands & _CURRENCY_NAMES:
        return True
    # Generic electron shuttles whose names vary between databases.
    return any(p in raw for p in _CURRENCY_PATTERNS)


# --- generic cofactor placeholders -> the host's real carrier ------------------
# ModelSEED/KEGG write many oxidation steps against an ABSTRACT redox pair — e.g.
# "Reduced flavin"/"Flavin" (rather than a concrete FADH2/FAD) or "Acceptor"/
# "Reduced acceptor". Those placeholders carry no cross-reference to any real
# compound, so the Canonicalizer cannot unify them with the host's pool and they
# enter the model as brand-new metabolites appearing in exactly ONE reaction.
# A dead-end metabolite must have zero flux at steady state, so the entire designed
# route is silently forced to zero — which is why the nootkatol/nootkatone pathway
# was found but could not carry any production flux. Worse, the database cannot
# rescue itself: every ModelSEED reaction PRODUCING "Reduced flavin" is another
# monooxygenase CONSUMING it, so the generic pair can never turn over.
#
# Mapping each placeholder onto the host's equivalent native carrier is what a
# metabolic engineer does by hand. Values are host id bases in order of preference.
# Every substitution is reported to the user (see `cofactor_substitutions`), because
# it is an interpretation of the database, not a fact stated by it.
_GENERIC_COFACTOR_ALIASES: Dict[str, tuple] = {
    # ModelSEED's generic NAD(P)+/NAD(P)H stand-in, used by 363 reactions each.
    "nad p or nop": ("nad", "nadp"),
    "nadh p or nop": ("nadh", "nadph"),
    "nadp or nop": ("nadp", "nad"),
    "nadph or nop": ("nadph", "nadh"),
    "reduced flavin": ("fadh2", "fmnh2"),
    "flavin reduced": ("fadh2", "fmnh2"),
    "flavin": ("fad", "fmn"),
    "flavin oxidized": ("fad", "fmn"),
    "oxidized flavin": ("fad", "fmn"),
    "reduced acceptor": ("fadh2", "nadh", "pqh2", "q8h2"),
    "acceptor": ("fad", "nad", "pq", "q8"),
    "oxidized acceptor": ("fad", "nad", "pq", "q8"),
    "reduced ferredoxin": ("fdxrd", "fdxr_42", "fdxr", "fdxred"),
    "oxidized ferredoxin": ("fdxox", "fdxo_42", "fdxo"),
    "ferredoxin": ("fdxox", "fdxo_42", "fdxo"),
    "reduced thioredoxin": ("trdrd", "trxrd"),
    "oxidized thioredoxin": ("trdox", "trxox"),
    "thioredoxin": ("trdox", "trxox"),
    "reduced plastocyanin": ("pcred", "pc_red"),
    "oxidized plastocyanin": ("pcox", "pc_ox"),
    "plastocyanin": ("pcox", "pc_ox"),
    "ferrocytochrome c": ("focytc", "focytC"),
    "ferricytochrome c": ("ficytc", "ficytC"),
    "cytochrome c reduced": ("focytc", "focytC"),
    "cytochrome c oxidized": ("ficytc", "ficytC"),
}


def _host_base_index(model: cobra.Model) -> Dict[str, cobra.Metabolite]:
    """Host metabolites indexed by lower-cased compartment-stripped id base."""
    idx: Dict[str, cobra.Metabolite] = {}
    for m in model.metabolites:
        idx.setdefault(_base_key(m.id).lower(), m)
    return idx


def resolve_generic_cofactor(met, host_index_by_base: Dict[str, cobra.Metabolite]):
    """Host carrier standing in for a generic cofactor placeholder, or ``None``.

    Returns ``(host_metabolite, generic_name)`` so the caller can report the swap.
    """
    raw = (getattr(met, "name", "") or "").strip().lower()
    if not raw:
        return None
    # Acyl thioesters are real carbon skeletons, never a cofactor placeholder.
    if ("coa" in raw and raw not in ("coa", "coenzyme a")) or "acp" in raw:
        return None
    # Singular and plural both, because databases disagree: the merged universal writes
    # "Reduced-ferredoxins" where the alias table (and BiGG) say "reduced ferredoxin".
    # Without this the lookup missed by one letter and the carrier entered the model as
    # a dead end, taking the whole route's flux to zero.
    candidates = _currency_name_candidates(raw) | {raw}
    candidates |= {c[:-1] for c in list(candidates) if c.endswith("s") and len(c) > 3}
    for cand in candidates:
        bases = _GENERIC_COFACTOR_ALIASES.get(cand)
        if not bases:
            continue
        for b in bases:
            hm = host_index_by_base.get(b.lower())
            if hm is not None:
                return hm, cand

    # The alias table names a carrier family by the ids it usually wears; a host that
    # spells one differently (iJN678's oxidised ferredoxin is `fdxo_2_2`, not `fdxo`)
    # would otherwise need a new literal added for every model. Fall back to matching
    # the cofactor *class* and redox state, which is model-independent.
    from . import id_reconcile
    wanted_class = id_reconcile.cofactor_class(met)
    if wanted_class:
        wanted_state = id_reconcile._redox_state(met)
        for hm in host_index_by_base.values():
            if id_reconcile.cofactor_class(hm) != wanted_class:
                continue
            # Never collapse the two halves of a couple into one another.
            if wanted_state and id_reconcile._redox_state(hm) != wanted_state:
                continue
            return hm, f"{wanted_state or ''} {wanted_class}".strip()
    return None


def currency_keys(*models) -> set:
    """Compartment-independent base keys of every currency/cofactor metabolite across
    the given models, detected namespace-agnostically (id AND name).

    The reachability search works on base-key strings, so it needs this set to know
    which substrates are freely available. Building it from the actual metabolite
    objects is what lets ModelSEED's "Reduced_flavin_c" / KEGG's "C00004" be
    recognised as cofactors rather than mistaken for carbon skeletons that would have
    to be synthesised first."""
    keys: set = set()
    for mdl in models:
        if mdl is None:
            continue
        for m in mdl.metabolites:
            if _is_currency_met(m):
                keys.add(_base_key(m.id).lower())
    return keys


def _carbon_count(met) -> int:
    """Number of carbon atoms in a metabolite's formula (0 if unknown). Used to
    pick the principal carbon skeleton in a pathway drawing."""
    formula = getattr(met, "formula", None) or ""
    total = 0
    for elem, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if elem == "C":
            total += int(num) if num else 1
    return total


def readable_metabolite_id(met: cobra.Metabolite) -> str:
    """A short human-readable id for a database metabolite (e.g. ``ibcoa_c``):
    from a BiGG id, else a slug of its name, else SEED/KEGG id, else its id."""
    import re
    comp = met.compartment or (met.id.rsplit("_", 1)[-1] if "_" in met.id else "c")
    ann = getattr(met, "annotation", None) or {}

    def _first(key):
        v = ann.get(key)
        if not v:
            return ""
        return (v[0] if isinstance(v, (list, tuple)) else str(v)).strip()

    bigg = _first("bigg.metabolite")
    if bigg:
        return f"{bigg}_{comp}"
    if met.name and met.name != met.id:
        slug = re.sub(r"[^A-Za-z0-9]+", "_", met.name).strip("_")[:24]
        if slug and not slug[0].isdigit():
            return f"{slug}_{comp}"
    for key in ("seed.compound", "kegg.compound"):
        val = _first(key)
        if val:
            return f"{val}_{comp}"
    return met.id


def _equation_with_names(rxn: cobra.Reaction) -> str:
    """Reaction equation using short/readable metabolite names."""
    from .network_graph import reaction_equation
    return reaction_equation(rxn)


def _find_pathway_reactions(host: cobra.Model, universal: cobra.Model, target_id: str,
                            max_steps: int, available_ids: Optional[set] = None,
                            preferred_ec: Optional[set] = None,
                            forbidden: Optional[set] = None,
                            require_balanced: bool = True,
                            include_boundary: bool = False) -> Optional[List[str]]:
    """Find a route to ``target_id`` by *network expansion* (scope).

    Starting from ``available_ids`` (default: all host metabolites) plus currency
    metabolites, reactions whose substrates are all reachable are "fired", adding
    their products — until the target is reachable; then the producing reactions
    are recovered by back-tracking. ``preferred_ec`` biases route choice toward
    reactions with those EC classes; ``forbidden`` reactions are skipped (used to
    enumerate alternative pathways).
    """
    from .databases import reaction_ec_numbers
    forbidden = forbidden or set()
    # Work in a compartment-independent key space so a reaction counts as the same
    # chemistry regardless of the compartment the database used, and pure transport
    # reactions (same species moved between compartments) contribute nothing.
    if available_ids is not None:
        avail = {_base_key(x) for x in available_ids}
    else:
        avail = {_base_key(m.id) for m in host.metabolites}
    target_key = _base_key(target_id)
    if target_key in avail:
        return []

    preferred_ec = preferred_ec or set()
    # Currency detection must be namespace-agnostic: the search works on id strings,
    # but only the metabolite *objects* carry the names that identify cofactors in
    # ModelSEED/MetaNetX (e.g. "Reduced_flavin_c", "cpd11621"). Resolve them once here
    # and test membership below, instead of pattern-matching BiGG-style ids.
    cur_keys = currency_keys(universal, host)

    def _cur(k: str) -> bool:
        return k.lower() in cur_keys or _is_currency_key(k)

    directions = []  # (pref, rxn_id, substrate_keys, product_keys, substrates_nc_keys)
    for r in universal.reactions:
        if r.id in forbidden:
            continue
        # By default, don't route through the database's own boundary reactions
        # (exchange/demand/sink) — they would let a metabolite appear "from nothing".
        if not include_boundary and r.boundary:
            continue
        # Issue 2: never route through a verifiably mass/charge-unbalanced reaction
        # (e.g. a C3 → C20 shortcut). Unverifiable reactions (missing formulas) are
        # kept but flagged downstream.
        if require_balanced and _is_grossly_unbalanced(r):
            continue
        subs = [_base_key(m.id) for m, c in r.metabolites.items() if c < 0]
        prods = [_base_key(m.id) for m, c in r.metabolites.items() if c > 0]
        subs_nc = {k for k in subs if not _cur(k)}
        prods_nc = {k for k in prods if not _cur(k)}
        if not prods_nc and not subs_nc:
            continue
        pref = 0
        if preferred_ec and (set(reaction_ec_numbers(r)) & preferred_ec):
            pref = 1
        directions.append((pref, r.id, subs, prods, subs_nc))
        if r.reversibility or (r.lower_bound < 0 < r.upper_bound):
            directions.append((pref, r.id, prods, subs, prods_nc))
    # try preferred-EC reactions first so routes favour them
    directions.sort(key=lambda d: -d[0])

    reachable = set(avail)
    pred: Dict[str, tuple] = {}
    found = target_key in reachable
    for _ in range(max(2, max_steps) + 4):
        changed = False
        for _pref, rid, subs, prods, subs_nc in directions:
            if subs_nc <= reachable:
                for p in prods:
                    if _cur(p) or p in reachable:
                        continue
                    reachable.add(p)
                    pred[p] = (rid, subs)
                    changed = True
        if target_key in reachable:
            found = True
            break
        if not changed:
            break
    if not found:
        return None

    ordered: List[str] = []
    seen_rxn = set()
    stack = [target_key]
    visited = set()
    while stack:
        k = stack.pop()
        if k in avail or k in visited or k not in pred:
            continue
        visited.add(k)
        rid, subs = pred[k]
        if rid not in seen_rxn:
            seen_rxn.add(rid)
            ordered.append(rid)
        for s in subs:
            if not _is_currency_key(s) and s not in avail:
                stack.append(s)
    ordered.reverse()
    return ordered


# How many extra routes to probe when every route found so far is infeasible
# (found, but unable to carry flux). Each probe costs one FBA, so keep it small.
_FLUX_RETRY_BUDGET = 4


def _carries_flux(r: "PathwayResult") -> bool:
    """Whether a result predicts real production (not NaN/infeasible, not zero)."""
    f = r.production_flux
    return f == f and f > 1e-9


def _n_carbons(met: cobra.Metabolite) -> int:
    """Carbon atoms in a metabolite: from its formula, else from its structure.

    The formula is missing far more often than one would expect — BiGG's universal model
    carries none at all — so fall back to any structure the metabolite is annotated with
    (the bundled database now ships SMILES/InChI where they could be resolved). Returns
    0 when neither is available, which callers must read as *unknown*, not as *zero*.
    Use :func:`_carbon_count` when that distinction has to survive.
    """
    return _carbon_count(met)[0]


def _carbon_count(met: cobra.Metabolite) -> tuple:
    """``(carbons, known)`` for a metabolite.

    ``known`` is the whole point: 0 carbons and "we have no idea" are different facts,
    and reading the second as the first is what turns a carbon yield into a silent NaN
    (or, worse, into a plausible-looking number computed over half the carbon).
    A metabolite with a real formula that contains no carbon is ``(0, True)``.
    """
    try:
        elements = met.elements or {}
    except Exception:  # noqa: BLE001 - malformed formula (R groups, polymers)
        elements = None
    if elements:
        n = int(elements.get("C", 0))
        if n:
            return n, True
    from_structure = _carbons_from_structure(met)
    if from_structure:
        return from_structure, True
    # No carbon found. Only trust that as a real zero when a formula actually parsed.
    parsed_a_formula = bool((getattr(met, "formula", "") or "").strip()) and elements
    return 0, bool(parsed_a_formula)


def _carbons_from_structure(met: cobra.Metabolite) -> int:
    ann = getattr(met, "annotation", None) or {}
    for key in ("SMILES", "smiles"):
        val = ann.get(key)
        if val:
            smiles = val[0] if isinstance(val, list) else val
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is not None:
                    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
            except Exception:  # noqa: BLE001 - RDKit absent or structure unparseable
                return 0
    for key in ("InChI", "inchi"):
        val = ann.get(key)
        if val:
            inchi = val[0] if isinstance(val, list) else val
            # The formula sits in the first layer: InChI=1S/C4H10O/c...
            parts = str(inchi).split("/")
            if len(parts) > 1:
                m = re.match(r"^C(\d*)(?![a-z])", parts[1])
                if m:
                    return int(m.group(1) or 1)
    return 0


@dataclass
class CarbonYield:
    """A carbon yield, or a plain-language reason there isn't one.

    NaN is not an answer. A yield that came out NaN because the product has no formula
    is a *database gap the user can fix*, and looks nothing like a yield that came out
    NaN because the route carries no flux — but as a bare float the two are identical.
    """

    value: float = float("nan")
    reason: str = ""            # empty exactly when `value` is a real number

    @property
    def computable(self) -> bool:
        return self.value == self.value          # not NaN

    def text(self) -> str:
        if self.computable:
            return f"{self.value * 100:.1f}% of consumed C"
        return f"not computable — {self.reason}" if self.reason else "not computable"


def _carbon_yield(model: cobra.Model, target_id: str, product_flux: float) -> CarbonYield:
    """mol C in the product / mol C taken up from the medium, at the current solution.

    Raw production flux is not comparable between products of different size, so this is
    the figure to use when ranking targets against each other (L10). Returns a
    :class:`CarbonYield` that says *why* when there is no number: a missing formula (the
    common case — roughly 40% of the merged database has none), no carbon exchange
    resolved, or an infeasible solution.
    """
    from . import balancing

    if not model.metabolites.has_id(target_id):
        return CarbonYield(reason=f"{target_id} is not in the engineered model")
    target = model.metabolites.get_by_id(target_id)
    prod_c, prod_known = _carbon_count(target)
    if not prod_known:
        return CarbonYield(reason=f"{balancing.metabolite_label(target)} has no formula")
    if prod_c <= 0:
        return CarbonYield(reason=f"{balancing.metabolite_label(target)} contains no carbon")
    if not (product_flux == product_flux) or product_flux <= 0:
        return CarbonYield(reason="the route carries no flux")
    sol = model.optimize()
    if sol.status != "optimal":
        return CarbonYield(reason="the engineered model has no optimal solution")
    uptake_c = 0.0
    unknown_uptake = []
    for rxn in model.reactions:
        if not rxn.boundary:
            continue
        v = float(sol.fluxes.get(rxn.id, 0.0))
        if v >= -1e-9:                     # only consumption (negative flux) is uptake
            continue
        for met in rxn.metabolites:
            c, known = _carbon_count(met)
            if not known:
                # Counting this uptake as carbon-free would inflate the yield — the same
                # absent-is-zero error, just hidden in the denominator this time.
                label = balancing.metabolite_label(met)
                if label not in unknown_uptake:
                    unknown_uptake.append(label)
            elif c:
                uptake_c += abs(v) * c
    if unknown_uptake:
        shown = ", ".join(unknown_uptake[:3])
        extra = len(unknown_uptake) - 3
        if extra > 0:
            shown += f" and {extra} other(s)"
        return CarbonYield(reason=f"consumed metabolite(s) with no formula: {shown}")
    if uptake_c <= 1e-9:
        return CarbonYield(reason="no carbon uptake was resolved in this solution")
    return CarbonYield(value=(product_flux * prod_c) / uptake_c)


def _result_from_ids(host, universal, target_id, ids, *, rank=1) -> PathwayResult:
    from .databases import reaction_ec_numbers
    rows = []
    warnings: List[str] = []
    unverified: List[str] = []       # steps whose balance could not be checked at all
    unverified_reasons: Dict[str, List[str]] = {}
    all_balanced = True
    for rid in ids:
        rxn = universal.reactions.get_by_id(rid)
        verdict = reaction_balance_verdict(rxn)
        if not verdict.checkable:
            # A missing formula is a DATABASE gap, not a chemical problem. Treating it as
            # an imbalance put a scary warning on almost every route, which trained users
            # to ignore the warning entirely — dangerous, because a real imbalance looks
            # identical. Unverifiable steps are recorded separately and do not raise the
            # headline warning (VI.14a). Name the culprit so the gap can be closed.
            bal_label = "unverified"
            unverified.append(rid)
            unverified_reasons[rid] = list(verdict.missing_formula_labels)
        elif verdict.verifiably_imbalanced:
            bal_label = f"UNBALANCED ({_format_residual(verdict.residual)})"
            warnings.append(f"{rid}: {verdict.sentence()}.")
            all_balanced = False
        else:
            bal_label = "balanced"
        # NOTE (#B4): the balance verdict is NOT shown as a table column any more —
        # it drives the result warning + the "Balance H+/H2O" button instead (bal_label
        # is retained only for the per-reaction right-click balancer).
        rows.append({
            "reaction": rxn.id,
            "suggested_id": readable_reaction_id(rxn),
            "name": rxn.name if rxn.name and rxn.name != rxn.id else "",
            "EC": ", ".join(reaction_ec_numbers(rxn)),
            "equation": _equation_with_names(rxn),
        })
    # Issue 2: only report a production flux when the engineered model can actually
    # carry it (FBA feasibility) — a physically impossible route yields NaN/0.
    flux = float("nan")
    apply_report: Dict[str, list] = {}
    engineered = None
    try:
        engineered = apply_pathway(host, universal, ids, report=apply_report,
                                   target=target_id)
        m = _ensure_metabolite(engineered, universal, target_id)
        demand_id = f"DM_{m.id}"
        dm = (engineered.reactions.get_by_id(demand_id) if engineered.reactions.has_id(demand_id)
              else engineered.add_boundary(m, type="demand"))
        engineered.objective = dm
        val = engineered.slim_optimize()
        flux = float(val) if val is not None and val == val else float("nan")
    except Exception:  # noqa: BLE001
        flux = float("nan")
    # A normalised, cross-target-comparable figure (L10). Raw flux cannot be compared
    # between products — a C2 and a C40 molecule make very different demands on the same
    # carbon budget — so also report carbon yield: C atoms leaving in the product divided
    # by C atoms entering from the medium at that solution.
    cy = CarbonYield(reason="the route carries no flux")
    try:
        if engineered is None:
            cy = CarbonYield(reason="the pathway could not be applied to the host")
        elif flux == flux and flux > 1e-9:
            cy = _carbon_yield(engineered, target_id, flux)
    except Exception:  # noqa: BLE001
        cy = CarbonYield(reason="the carbon balance could not be evaluated")
    prefix = f"Pathway {rank}: " if rank else ""
    note = f"{prefix}{len(ids)} heterologous reaction(s) to {target_id}."
    if not all_balanced:
        # Reserved for genuine mass/charge imbalance — see VI.14a.
        note += (" ⚠ Some step(s) are mass/charge-unbalanced. Check them before relying "
                 "on this route (the “Balance H⁺/H₂O” button fixes the common "
                 "protonation/hydration cases).")
    # A route can be perfectly good chemistry yet carry no flux, because the search is
    # topological and does not check that by-products can be disposed of. Say WHY, so
    # the user sees an actionable cause instead of an unexplained zero.
    swaps = apply_report.get("cofactor_substitutions") or []
    if swaps:
        warnings.append(
            "The database left some cofactors generic; your model's carrier was used "
            "instead: " + ", ".join(swaps) + ". Check these are the intended cofactors.")
    dead = apply_report.get("dead_end_metabolites") or []
    blocked_by: List[str] = []
    if dead and not (flux == flux and flux > 1e-9):
        shown = ", ".join(dead[:6]) + (" …" if len(dead) > 6 else "")
        note += (" This route cannot carry flux as it stands: " + shown + " can only be "
                 "produced or only consumed, so at steady state the pathway is blocked. "
                 "Add a reaction consuming/producing it (or an exchange), or try another "
                 "route.")
        warnings.append(f"Dead-end metabolite(s) block this route: {shown}")
        # Pin down WHICH step introduced the dead end, so the search can be retried
        # without just that step. `apply_pathway` renames reactions, so map the
        # engineered ids back to the database ids they came from.
        if engineered is not None:
            eng_to_db = {}
            for rid in ids:
                src = universal.reactions.get_by_id(rid)
                eng_to_db[src.id] = rid
                eng_to_db[readable_reaction_id(src)] = rid
            for mid in dead:
                if not engineered.metabolites.has_id(mid):
                    continue
                for r in engineered.metabolites.get_by_id(mid).reactions:
                    db_rid = eng_to_db.get(r.id)
                    if db_rid and db_rid not in blocked_by:
                        blocked_by.append(db_rid)
    # Skeletal-isomer check (VI.1A). Offline only here — it must not slow a search down
    # with network calls; the GUI re-runs it with `online=True` for the feasibility report.
    iso_warnings: List[str] = []
    iso_checked = False
    iso_coverage = ""
    try:
        from .chemistry import check_pathway_isomers
        rep = check_pathway_isomers(universal, ids, online=False)
        iso_checked = bool(rep.checked)
        iso_coverage = rep.coverage
        iso_warnings = [f.as_sentence() for f in rep.findings]
        if iso_warnings:
            warnings.extend(iso_warnings)
            note += (" 🚩 A step changes the carbon backbone in a way no single reaction "
                     "can — this route may be mislabelled chemistry. See Feasibility "
                     "information.")
    except Exception:  # noqa: BLE001 — a check must never break a search
        pass

    return PathwayResult(
        target=target_id, reactions=pd.DataFrame(rows), n_steps=len(ids),
        production_flux=flux, reaction_ids=ids, note=note, blocked_by=blocked_by,
        warnings=warnings, balanced=all_balanced, carbon_yield=cy.value,
        carbon_yield_note=cy.reason,
        unverified_steps=unverified, unverified_reasons=unverified_reasons,
        isomer_warnings=iso_warnings,
        isomer_checked=iso_checked, isomer_coverage=iso_coverage)


# Redox/energy carriers whose identity DOES distinguish two otherwise identical routes:
# an NADH-dependent and an NADPH-dependent variant are a real engineering choice, so they
# must stay separate even though the carbon chemistry is the same (VI.12).
_COFACTOR_DISCRIMINATORS = (
    ("nadh", "nad_"), ("nadph", "nadp_"), ("fadh2", "fad_"), ("fmnh2", "fmn_"),
    ("atp", "adp"), ("q8h2", "q8"), ("mql8", "mqn8"),
    ("ferredoxin", "flavodoxin"), ("cytochrome",),
)


def _reaction_chemistry_key(rxn) -> tuple:
    """What a reaction DOES, ignoring which database entry it came from.

    Two entries that convert the same carbon compounds are the same chemistry even when
    their ids differ (``rxn35741`` vs ``alcohol_dehydrogenase_NA_31``). Currency
    metabolites are excluded from the key so the comparison is about the substrate and
    product, not the cofactors.
    """
    from .network_graph import CURRENCY_BASES
    subs, prods = [], []
    for met, coeff in rxn.metabolites.items():
        base = met.id.rsplit("_", 1)[0].lower() if "_" in met.id else met.id.lower()
        if base in CURRENCY_BASES:
            continue
        (subs if coeff < 0 else prods).append(base)
    return (tuple(sorted(subs)), tuple(sorted(prods)))


def _cofactor_signature(rxn) -> tuple:
    """Which redox/energy carriers a reaction uses — the part that must NOT be collapsed.

    Matching is on the metabolite's *base id* (compartment stripped), not a substring,
    so ``nadp_c`` is not mistaken for ``adp``.
    """
    bases = set()
    for m in rxn.metabolites:
        mid = m.id.lower()
        bases.add(mid.rsplit("_", 1)[0] if "_" in mid else mid)
    sig = []
    for group in _COFACTOR_DISCRIMINATORS:
        for token in group:
            tok = token.rstrip("_")
            # exact base match, or a prefix match for the named families (cytochrome…)
            if tok in bases or any(b == tok or b.startswith(tok) and len(tok) > 5
                                   for b in bases):
                sig.append(tok)
    return tuple(sorted(set(sig)))


def _route_signature(result, universal) -> tuple:
    """A route's identity: its chemistry plus the cofactors it commits you to."""
    chem, cof = [], []
    for rid in result.reaction_ids:
        if not universal.reactions.has_id(rid):
            chem.append((rid,))
            continue
        rxn = universal.reactions.get_by_id(rid)
        chem.append(_reaction_chemistry_key(rxn))
        cof.append(_cofactor_signature(rxn))
    return (tuple(sorted(map(str, chem))), tuple(sorted(map(str, cof))))


def _collapse_equivalent_routes(results: List[PathwayResult], universal
                                ) -> List[PathwayResult]:
    """Merge alternatives that are the same chemistry with the same cofactors (VI.12).

    The alternatives counter used to overstate how much genuine choice the user had:
    three of six "different" n-butanol routes were one reaction under three database ids.
    Routes differing only in NADH vs NADPH are kept apart, because that is a real design
    decision, and the surviving route records the ids it stands for.
    """
    seen: Dict[tuple, PathwayResult] = {}
    out: List[PathwayResult] = []
    for r in results:
        if not r.reaction_ids:
            out.append(r)
            continue
        try:
            sig = _route_signature(r, universal)
        except Exception:  # noqa: BLE001 — never lose a result to a signature failure
            out.append(r)
            continue
        keep = seen.get(sig)
        if keep is None:
            seen[sig] = r
            out.append(r)
            continue
        # Same chemistry AND same cofactors: keep the better one, remember the other.
        better = r if (_carries_flux(r) and not _carries_flux(keep)) else keep
        other = keep if better is r else r
        better.equivalent_routes = list(
            getattr(better, "equivalent_routes", []) or []) + list(other.reaction_ids)
        if better is r:                       # swap the kept route in place
            out[out.index(keep)] = r
            seen[sig] = r
    for r in out:
        eq = getattr(r, "equivalent_routes", None)
        if eq:
            r.note += (f" ({len(set(eq))} equivalent database entr"
                       f"{'y' if len(set(eq)) == 1 else 'ies'} collapsed — same chemistry "
                       "and cofactors.)")
    return out


def find_pathways(
    host: cobra.Model,
    target_metabolite_id: str,
    universal: cobra.Model,
    *,
    start_metabolites: Optional[List[str]] = None,
    max_steps: int = 25,
    n_alternatives: int = 1,
    preferred_ec: Optional[List[str]] = None,
    forbidden_reactions: Optional[set] = None,
    match_namespace: bool = True,
    algorithm: str = "retro",
    strategy: str = "expansion",           # deprecated alias; kept for callers
    priority: str = "yield",
    include_boundary: bool = False,
    progress: Optional[Callable[[int, int, str], None]] = None,
    collapse_duplicates: bool = True,
) -> List[PathwayResult]:
    """Return up to ``n_alternatives`` distinct heterologous routes to the target,
    using the rebuilt identity-aware search engine (:mod:`pathway_search`).

    ``algorithm`` selects the search:

    * ``"retro"`` (default) — retrosynthetic AND/OR search that finds the *minimal*
      set of heterologous reactions connecting the target back to compounds the
      host already makes (native precursors end a branch immediately).
    * ``"expansion"`` — forward network expansion from the host's compounds.

    Both are compartment-agnostic and unify host/database compounds by their
    cross-references, so a native precursor with a different database id is
    recognised. Alternatives are found by forbidding a route's reactions and
    re-searching; ``priority`` (``"yield"``/``"shortest"``) ranks them.
    """
    from . import pathway_search

    pref_ec = set(preferred_ec or [])
    # NOTE: `start_metabolites` must NOT imply `include_boundary`. Naming the allowed
    # precursors is a *narrowing* constraint, but ORing it into include_boundary used to
    # switch ON the database's own exchange/demand reactions — whose substrate set is
    # empty, so they fire unconditionally and let ANY compound with an exchange appear
    # from nothing. That is what made routes rely on an external supply of a compound
    # the host cannot make (e.g. hexanoyl-CoA in Synechocystis). Boundary reactions are
    # now used only when the user explicitly ticks the option, and the start list is
    # passed to the search to restrict the precursor pool as intended.

    # Target must exist somewhere (host or database), in any compartment.
    if pathway_search._resolve_target(host, universal, target_metabolite_id) is None:
        return [PathwayResult(
            target=target_metabolite_id, production_flux=0.0,
            note=f"'{target_metabolite_id}' is not in the model or the reaction database.")]

    def _search(forbidden):
        if algorithm == "expansion":
            return pathway_search.expansion_search(
                host, universal, target_metabolite_id, max_steps=max_steps,
                preferred_ec=pref_ec, forbidden=forbidden, include_boundary=include_boundary,
                start_metabolites=start_metabolites)
        return pathway_search.retro_search(
            host, universal, target_metabolite_id, max_reactions=max_steps,
            preferred_ec=pref_ec, forbidden=forbidden, include_boundary=include_boundary,
            start_metabolites=start_metabolites)

    results: List[PathwayResult] = []
    forbidden: set = set(forbidden_reactions or set())
    wanted = max(1, n_alternatives)
    if progress:
        progress(0, wanted, "Searching for pathway 1…")
    # The search is topological: it guarantees every substrate is reachable, but not
    # that the route's by-products can be disposed of. Such a route is real chemistry
    # yet carries ZERO flux at steady state, which is useless to the user. So if the
    # routes found so far are all infeasible, spend a small extra budget looking for
    # one that actually runs instead of returning the first dead one.
    for _try in range(wanted + _FLUX_RETRY_BUDGET):
        ids = _search(forbidden)
        if ids is None:
            break
        if not ids:  # host already makes the target
            results.append(PathwayResult(
                target=target_metabolite_id, production_flux=float("nan"),
                note="The host already produces this metabolite — no heterologous "
                     "reactions are needed."))
            break
        res = _result_from_ids(host, universal, target_metabolite_id, ids, rank=0)
        results.append(res)
        if progress:
            # Report against the number the user asked for, not the retry budget, so the
            # dialog reads "pathway 2/4" rather than a confusing larger total (VI.14b).
            done = min(len(results), wanted)
            progress(done, wanted,
                     f"Found pathway {done} of {wanted}…" if done < wanted
                     else "Finishing…")
        if _carries_flux(res) or not res.blocked_by:
            forbidden |= set(ids)         # enumerate a genuinely different alternative
        else:
            # This route is blocked by a dead end. Forbid ONLY the step that causes it:
            # banning the whole route would also ban the steps it shares with the
            # working alternative — including the essential terminal enzyme — leaving
            # no route findable at all.
            forbidden |= set(res.blocked_by)
        if len(results) >= wanted and any(_carries_flux(r) for r in results):
            break

    if collapse_duplicates and len(results) > 1:
        results = _collapse_equivalent_routes(results, universal)

    # Issue 3/R6: rank alternatives by the chosen priority — highest predicted
    # production flux ("yield") or fewest heterologous reactions ("shortest").
    # Balanced routes rank above unbalanced ones; NaN (infeasible) fluxes sink.
    # A route that carries flux ALWAYS outranks one that cannot run at all, whatever
    # the priority — a shorter dead route is not a better answer. Zero-flux routes are
    # RANKED DOWN, never discarded: they are still real chemistry and often the only
    # option, so the user should see them (with the reason attached) and decide.
    if len(results) > 1 and any(r.reaction_ids for r in results):
        def _key(r):
            f = r.production_flux
            f = f if (f == f) else float("-inf")  # NaN -> last
            if priority == "shortest":
                return (_carries_flux(r), r.balanced, -r.n_steps, f)
            return (_carries_flux(r), r.balanced, f, -r.n_steps)
        results.sort(key=_key, reverse=True)
        del results[max(1, n_alternatives):]      # trim the extra feasibility probes
        label = "fewest steps" if priority == "shortest" else "yield"
        for i, r in enumerate(results, start=1):
            if r.reaction_ids:
                r.note = f"Pathway {i} (ranked by {label}): " + r.note.split(": ", 1)[-1]

    if not results:
        # "0 pathways" is not an answer. It conflates three different situations — the
        # compound is absent; it is catalogued but nothing produces it; or a precursor
        # is unreachable — and only the last is helped by allowing more steps. Say
        # which, and name the compound that breaks the chain.
        note = f"No heterologous pathway to {target_metabolite_id} was found."
        warnings: List[str] = []
        missing: List[str] = []
        try:
            from . import pathway_diagnostics
            why = pathway_diagnostics.explain_not_found(host, universal,
                                                        target_metabolite_id)
            if why.summary:
                note = why.summary
            if why.recommendation:
                note += " " + why.recommendation
            if why.missing:
                missing = list(why.missing)
                warnings.append("No reaction in this database produces: "
                                + ", ".join(why.missing[:4]))
        except Exception:  # noqa: BLE001 — a diagnosis must never break the search
            shared = sum(1 for m in universal.metabolites if host.metabolites.has_id(m.id))
            note += (" The database shares no metabolites with your model — try one in "
                     "your model's namespace." if shared == 0 else
                     " Try increasing the step limit or including more databases.")
        return [PathwayResult(target=target_metabolite_id, production_flux=0.0,
                              note=note, warnings=warnings,
                              missing_compounds=missing)]
    return results


def predict_pathway(
    host: cobra.Model,
    target_metabolite_id: str,
    universal: cobra.Model,
    *,
    lower_bound: float = 0.1,
    iterations: int = 1,
    match_namespace: bool = True,
    max_steps: int = 12,
) -> PathwayResult:
    """Predict a single heterologous route to a target (see :func:`find_pathways`)."""
    return find_pathways(host, target_metabolite_id, universal, max_steps=max_steps,
                         n_alternatives=1, match_namespace=match_namespace)[0]


def readable_reaction_id(rxn: cobra.Reaction) -> str:
    """A short, human-friendly id for a database reaction (e.g. ``PAL``).

    Preference: a genuine short id (BiGG/``short_id``), then a slug of the
    descriptive reaction *name* (so opaque accessions like ``R12101``/``MNXR…`` are
    a last resort), then those accessions, then the raw id."""
    ann = getattr(rxn, "annotation", None) or {}

    def _first(key):
        val = ann.get(key)
        if not val:
            return ""
        return (val[0] if isinstance(val, (list, tuple)) else str(val)).strip()

    def _name_slug():
        if rxn.name and rxn.name != rxn.id:
            slug = re.sub(r"[^A-Za-z0-9]+", "_", rxn.name).strip("_")[:24]
            if slug and not slug[0].isdigit():
                return slug
        return ""

    for key in ("short_id", "bigg.reaction"):
        v = _first(key)
        if v:
            return v
    slug = _name_slug()
    if slug:
        return slug
    for key in ("kegg.reaction", "seed.reaction"):
        v = _first(key)
        if v:
            return v
    return rxn.id


_RAW_MET_ID = re.compile(r"^(C\d{5}|MNXM\d+|cpd\d+|G\d{5}|CHEBI[:_]?\d+)(_[a-z0-9]+)?$", re.I)
_RAW_RXN_ID = re.compile(r"^(R\d{5}|MNXR\d+|rxn\d+)$", re.I)


def result_to_dict(r: "PathwayResult") -> dict:
    """Serialise a PathwayResult to a plain dict (for saving a designed pathway)."""
    flux = None if r.production_flux != r.production_flux else float(r.production_flux)
    reactions = r.reactions.to_dict("records") if not r.reactions.empty else []
    return {"target": r.target, "reaction_ids": list(r.reaction_ids),
            "n_steps": int(r.n_steps), "production_flux": flux, "note": r.note,
            "reactions": reactions}


def result_from_dict(d: dict) -> "PathwayResult":
    """Rebuild a PathwayResult from :func:`result_to_dict` output."""
    flux = float("nan") if d.get("production_flux") is None else float(d["production_flux"])
    return PathwayResult(
        target=d.get("target", ""), reactions=pd.DataFrame(d.get("reactions") or []),
        n_steps=int(d.get("n_steps", 0)), production_flux=flux,
        note=d.get("note", ""), reaction_ids=list(d.get("reaction_ids") or []))


def pathway_chain(universal: cobra.Model, reaction_ids: List[str], target_id: str):
    """Order a heterologous pathway as a linear chain for a scheme drawing.

    Returns ``(metabolite_ids, steps)`` where ``metabolite_ids`` is the sequence of
    *main* (non-currency) metabolites from the nearest precursor down to the target,
    and ``steps[i]`` is ``(reaction, consumed_currency_ids, produced_currency_ids)``
    converting ``metabolite_ids[i]`` → ``metabolite_ids[i+1]``. Reactions may be used
    in reverse (the currency lists are oriented accordingly). Falls back to the raw
    reaction order if the target cannot be back-tracked.
    """
    rxns = [universal.reactions.get_by_id(r) for r in reaction_ids
            if universal.reactions.has_id(r)]
    met_of = {m.id: m for rxn in rxns for m in rxn.metabolites}

    def _mains(rxn, produced_side: bool):
        # Main (carbon-skeleton) metabolites on one side, currency excluded by NAME
        # so it works for formula-less MetaNetX/KEGG metabolites.
        return [m for m, c in rxn.metabolites.items()
                if ((c > 0) == produced_side) and not _is_currency_met(m)]

    # Map each producible main metabolite -> list of (reaction, forward?) producing it.
    by_product: Dict[str, list] = {}
    for rxn in rxns:
        for m in _mains(rxn, True):
            by_product.setdefault(m.id, []).append((rxn, True))
        if rxn.reversibility:
            for m in _mains(rxn, False):
                by_product.setdefault(m.id, []).append((rxn, False))

    # If the target id isn't a product here (e.g. it was translated to a host id),
    # fall back to a terminal product: produced by the set but not consumed within it.
    if target_id not in by_product:
        consumed_all = set()
        for rxn in rxns:
            for m in _mains(rxn, False):
                consumed_all.add(m.id)
            if rxn.reversibility:
                for m in _mains(rxn, True):
                    consumed_all.add(m.id)
        sinks = [mid for mid in by_product if mid not in consumed_all]
        if sinks:
            target_id = sinks[0]

    # Back-track from the target to the precursor, using EACH reaction at most once
    # (a linear pathway uses each heterologous step once). At every step choose the
    # substrate that CONTINUES the chain — one made by a *different* unused reaction —
    # rather than a co-substrate that only exists by reversing the current step (which
    # is what made malonyl-CoA look like an intermediate).
    seq = []          # (rxn, forward, product_id, substrate_id)
    used: set = set()
    visited: set = set()
    cur = target_id
    while cur is not None and cur not in visited:
        visited.add(cur)
        producers = [(rxn, fwd) for (rxn, fwd) in by_product.get(cur, [])
                     if rxn.id not in used]
        if not producers:
            break
        rxn, forward = producers[0]
        used.add(rxn.id)
        subs = [m for m in _mains(rxn, not forward) if m.id not in visited] \
            or _mains(rxn, not forward)

        def _score(m):
            made_by_other = any(r.id not in used for (r, _f) in by_product.get(m.id, []))
            return (made_by_other, _carbon_count(m))

        nxt = max(subs, key=_score).id if subs else None
        seq.append((rxn, forward, cur, nxt))
        cur = nxt
    seq.reverse()

    if not seq:
        return [], []
    met_ids = [seq[0][3]] + [s[2] for s in seq]
    steps = []
    for rxn, forward, _prod, _sub in seq:
        cons_side, prod_side = (-1, 1) if forward else (1, -1)
        consumed = [m.id for m, c in rxn.metabolites.items()
                    if (c < 0) == (cons_side < 0) and _is_currency_met(m)]
        produced = [m.id for m, c in rxn.metabolites.items()
                    if (c < 0) == (prod_side < 0) and _is_currency_met(m)]
        # Extra main (carbon-skeleton) substrates feeding this step besides the chain
        # one — e.g. a lyase/synthase fusing two backbones. Drawn as merge inputs.
        cosubs = [m.id for m in _mains(rxn, not forward)
                  if _base_key(m.id) != _base_key(_sub or "")]
        steps.append((rxn, consumed, produced, cosubs))
    return met_ids, steps


def load_readable_cached(source_path: str, builder) -> cobra.Model:
    """Return a readable-id universal for ``source_path``, using a fast pickle cache.

    The first time, ``builder()`` produces the model (download/parse), it is made
    id-readable, and the finished model is pickled next to ``source_path``. On later
    runs the pickle is loaded directly — far faster than re-parsing a large JSON
    universal and re-running the id renaming every launch. The pickle is rebuilt if
    the source file is newer (or the pickle is missing/unreadable).
    """
    import pickle

    pkl = source_path + ".readable.pkl"
    try:
        if (os.path.exists(pkl) and os.path.exists(source_path)
                and os.path.getmtime(pkl) >= os.path.getmtime(source_path)):
            with open(pkl, "rb") as fh:
                return pickle.load(fh)
    except Exception:  # noqa: BLE001 - stale/corrupt cache -> rebuild
        pass
    model = builder()
    make_model_ids_readable(model)
    try:
        with open(pkl, "wb") as fh:
            pickle.dump(model, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception:  # noqa: BLE001 - caching is best-effort
        pass
    return model


def make_model_ids_readable(model: cobra.Model, *, reactions: bool = True) -> cobra.Model:
    """Rename raw database accessions (KEGG ``C08580``, MetaNetX ``MNXM…``) to
    human-readable ids in place, preserving the original id as an annotation.

    Ids that are already human-readable (e.g. BiGG ``glc__D_c``) are left alone,
    and renames never collide with an existing id. Cross-reference annotations are
    kept so the database still matches the host model across namespaces.

    Performance: the ids are assigned via the private ``_id`` attribute and the
    model index is rebuilt **once** at the end (``model.repair()``). Using the
    public ``.id`` setter on a model-attached object rebuilds the DictList index on
    *every* assignment — O(n²), which took ~6 minutes on a genome-scale MetaNetX
    universal; this makes it a few seconds.
    """
    used_met = {m.id for m in model.metabolites}
    renamed = False
    for met in model.metabolites:
        if not _RAW_MET_ID.match(met.id):
            continue
        new_id = readable_metabolite_id(met)
        if new_id == met.id or new_id in used_met:
            continue
        used_met.discard(met.id)
        used_met.add(new_id)
        met._id = new_id
        renamed = True
    if reactions:
        used_rxn = {r.id for r in model.reactions}
        for rxn in model.reactions:
            if not _RAW_RXN_ID.match(rxn.id):
                continue
            new_id = readable_reaction_id(rxn)
            # readable_reaction_id may return the raw accession itself (e.g. a KEGG
            # ``R…`` id); in that case prefer a slug of the reaction name.
            if _RAW_RXN_ID.match(new_id) and rxn.name and rxn.name != rxn.id:
                slug = re.sub(r"[^A-Za-z0-9]+", "_", rxn.name).strip("_")[:24]
                if slug:
                    new_id = slug
            if new_id == rxn.id or new_id in used_rxn:
                continue
            used_rxn.discard(rxn.id)
            used_rxn.add(new_id)
            rxn._id = new_id
            renamed = True
    if renamed:
        model.repair()
    return model


def add_pathway_to_model(model: cobra.Model, universal: cobra.Model,
                         reaction_ids: List[str], *, match_namespace: bool = False,
                         rename: bool = False, report: Optional[dict] = None,
                         rename_map: Optional[dict] = None,
                         target: Optional[str] = None, ensure_exchange: bool = True,
                         gas_product: bool = False,
                         product_route: str = "secrete") -> cobra.Model:
    """Add the given universal reactions (and any missing metabolites) to ``model`` in place.

    With ``match_namespace`` the universal is first translated into ``model``'s
    metabolite namespace, so the added reactions connect to the existing network
    even when the database used different identifiers. With ``rename`` the added
    reactions get a short human-readable id (e.g. ``PAL``) instead of the raw
    database id, when one is available and free.

    If a ``report`` dict is supplied it is filled (Issue 7) with keys
    ``added_reactions``, ``dropped_reactions`` (pure transport after collapse),
    ``collapsed_compartments`` and ``unbalanced`` (``[(rid, residual)]``) so the
    caller can warn the user about silent merges / imbalance after applying.
    """
    # Cross-namespace compound identity so a database metabolite that is really a
    # host-native compound (different id/compartment) is REUSED, not duplicated as
    # a new (green) metabolite — the root of the "natives turn green" bug (#B7).
    from . import pathway_search
    canon = pathway_search.Canonicalizer([model, universal])
    host_index: Dict[str, cobra.Metabolite] = {}
    for hm in model.metabolites:
        host_index.setdefault(canon.key(hm), hm)
    # Every metabolite the host had BEFORE this pathway was applied. `host_index` keeps
    # only one metabolite per canonical key, so it must not be used for this: the same
    # compound in another compartment (cu2_u vs cu2_c) would look newly added.
    native_ids = {m.id for m in model.metabolites}

    rep_added: List[str] = []
    rep_dropped: List[str] = []
    rep_collapsed: set = set()
    rep_cofactors: List[str] = []
    rep_dead_cofactors: List[str] = []
    # Index for mapping the database's generic redox placeholders ("Reduced flavin"…)
    # onto the host's real carrier, so they don't enter as dead-end metabolites that
    # force the whole route to zero flux.
    host_by_base = _host_base_index(model)

    # Compartment collapse: the host may lack compartments the database uses
    # (mitochondria, chloroplast…). We treat a reaction as its chemistry regardless
    # of compartment, mapping any metabolite in a compartment the host doesn't have
    # into the host's main compartment (usually the cytosol). Pure transport
    # reactions (a species moved between compartments) then cancel and are dropped.
    host_comps = set(getattr(model, "compartments", {}) or {})
    if not host_comps:
        host_comps = {m.compartment for m in model.metabolites if m.compartment}
    if "c" in host_comps:
        main_comp = "c"
    else:
        from collections import Counter
        cnt = Counter(m.compartment for m in model.metabolites if m.compartment)
        main_comp = cnt.most_common(1)[0][0] if cnt else "c"

    used_ids = {m.id for m in model.metabolites}
    met_map: Dict[str, cobra.Metabolite] = {}

    def _resolve(met: cobra.Metabolite) -> cobra.Metabolite:
        if met.id in met_map:
            return met_map[met.id]
        # 1) Is this compound already in the host (by identity)? Reuse it — this is
        #    a native metabolite and must NOT become a new/green one.
        hm = host_index.get(canon.key(met))
        if hm is not None:
            met_map[met.id] = hm
            return hm
        # 1b) A generic cofactor placeholder the database never resolves (e.g.
        #     "Reduced flavin"/"Flavin"): stand the host's real carrier in for it.
        #     Otherwise it becomes a dead-end compound and the route carries no flux.
        alias = resolve_generic_cofactor(met, host_by_base)
        if alias is not None:
            hm, generic = alias
            met_map[met.id] = hm
            rep_cofactors.append(f"{met.name or met.id} → {hm.id}")
            return hm
        comp = met.compartment or (met.id.rsplit("_", 1)[-1] if "_" in met.id else main_comp)
        new_comp = comp if comp in host_comps else main_comp
        if new_comp != comp:
            rep_collapsed.add(f"{comp}→{new_comp}")
        readable = readable_metabolite_id(met) if rename else met.id
        base = readable.rsplit("_", 1)[0] if "_" in readable else readable
        cand = f"{base}_{new_comp}"
        # 2) Reuse an existing host metabolite by id, if present.
        for cid in (met.id, cand):
            if model.metabolites.has_id(cid):
                m = model.metabolites.get_by_id(cid)
                met_map[met.id] = m
                return m
        # 3) Genuinely new (heterologous) compound — create it, cytosolic.
        new_id = cand if cand not in used_ids else met.id
        m = cobra.Metabolite(new_id, name=met.name or new_id, formula=met.formula,
                             charge=met.charge, compartment=new_comp)
        m.annotation = dict(getattr(met, "annotation", {}) or {})
        m.annotation.setdefault("metanetx.chemical", met.id.rsplit("_", 1)[0])
        used_ids.add(new_id)
        met_map[met.id] = m
        return m

    to_add = []
    for rid in reaction_ids:
        if not universal.reactions.has_id(rid):
            continue
        src = universal.reactions.get_by_id(rid)
        # A user-edited "suggested id" (#B5) wins over the auto-readable id.
        user_id = (rename_map or {}).get(src.id)
        new_id = user_id or (readable_reaction_id(src) if rename else src.id)
        if model.reactions.has_id(new_id) or model.reactions.has_id(src.id):
            continue
        coeffs: Dict[cobra.Metabolite, float] = {}
        for met, c in src.metabolites.items():
            rm = _resolve(met)
            coeffs[rm] = coeffs.get(rm, 0.0) + c
        coeffs = {m: c for m, c in coeffs.items() if abs(c) > 1e-9}
        if not coeffs:      # pure transport / identity after collapse -> not a real step
            rep_dropped.append(src.id)
            continue
        new = cobra.Reaction(new_id, name=src.name, subsystem=getattr(src, "subsystem", ""))
        new.bounds = src.bounds
        new.annotation = dict(getattr(src, "annotation", {}) or {})
        new.annotation.setdefault("metanetx.reaction", src.id)
        to_add.append((new, coeffs))
    model.add_reactions([nr for nr, _ in to_add])
    for new, coeffs in to_add:
        new.add_metabolites(coeffs)
        rep_added.append(new.id)


    # Streamline the expert workflow: guarantee the product can actually leave the
    # model (a single exchange/secretion route), so FBA/envelope/strain-design work
    # on the target straight away (#F3). Reuses an existing exchange if present.
    exchange_info = None
    if ensure_exchange and target:
        tmet = _resolve_target_metabolite(model, universal, target, met_map, host_index, canon)
        if tmet is not None:
            try:
                from .analysis.strain_design import ensure_product_exchange
                ex_id, created, kind = ensure_product_exchange(
                    model, tmet.id, route=product_route, gas=gas_product)
                exchange_info = {"id": ex_id, "created": bool(created),
                                 "metabolite": tmet.id, "route": kind}
            except Exception:  # noqa: BLE001 - never fail the apply over the exchange
                exchange_info = None

    # Safety net (runs AFTER the product exchange exists, or the target itself would
    # count as a dead end). The search is purely topological: it checks that every
    # substrate can be reached, but never that the route's by-products can be disposed
    # of, nor that a cofactor it waved through as "free" can actually be regenerated
    # here. Either gap leaves a metabolite that only ever appears on one side of the
    # reactions touching it — a dead end, which at steady state pins the whole route to
    # ZERO flux. Detect them so the user is told WHY a pathway cannot run.
    for m in list(model.metabolites):
        if m.id in native_ids:
            continue                       # native: the host's own network balances it
        makes = consumes = False
        for r in m.reactions:
            c = r.metabolites[m]
            fwd, rev = r.upper_bound > 0, r.lower_bound < 0
            if (c > 0 and fwd) or (c < 0 and rev):
                makes = True
            if (c < 0 and fwd) or (c > 0 and rev):
                consumes = True
            if makes and consumes:
                break
        if not (makes and consumes):
            rep_dead_cofactors.append(m.id)

    if report is not None:
        unbalanced = []
        unverified = []
        for rid in rep_added:
            rxn = model.reactions.get_by_id(rid)
            verdict = reaction_balance_verdict(rxn)
            if verdict.verifiably_imbalanced:
                unbalanced.append((rid, _format_residual(verdict.residual)))
            elif not verdict.checkable:
                # Not an imbalance — a gap in the database. Reported on its own so the
                # summary never accuses a reaction of losing mass it never counted.
                unverified.append((rid, verdict.missing_formula_sentence()))
        report["added_reactions"] = rep_added
        report["dropped_reactions"] = rep_dropped
        report["collapsed_compartments"] = sorted(rep_collapsed)
        report["unbalanced"] = unbalanced
        report["unverified"] = unverified
        report["product_exchange"] = exchange_info
        report["cofactor_substitutions"] = sorted(set(rep_cofactors))
        report["dead_end_metabolites"] = sorted(set(rep_dead_cofactors))
    return model


def _resolve_target_metabolite(model, universal, target, met_map, host_index, canon):
    """Find the host metabolite the target refers to after a pathway is applied,
    across id/namespace differences (used to attach the product exchange)."""
    if not target:
        return None
    if model.metabolites.has_id(target):
        return model.metabolites.get_by_id(target)
    if target in met_map:                       # universal id resolved during apply
        return met_map[target]
    if universal.metabolites.has_id(target):
        um = universal.metabolites.get_by_id(target)
        return met_map.get(um.id) or host_index.get(canon.key(um))
    base = target.rsplit("_", 1)[0] if "_" in target else target
    for m in model.metabolites:                 # last resort: base-id match
        if (m.id.rsplit("_", 1)[0] if "_" in m.id else m.id) == base:
            return m
    return None


def apply_pathway(host: cobra.Model, universal: cobra.Model, reaction_ids: List[str],
                  *, match_namespace: bool = False, rename: bool = False,
                  rename_map: Optional[dict] = None, report: Optional[dict] = None,
                  target: Optional[str] = None, ensure_exchange: bool = True,
                  gas_product: bool = False, product_route: str = "secrete") -> cobra.Model:
    """Return a copy of ``host`` with the given universal reactions added.

    ``rename_map`` (user-edited suggested ids), ``report`` (filled with the
    added/collapsed/dropped/unbalanced/product-exchange details), ``target``
    (the product metabolite, for which a single export route is guaranteed) and
    ``product_route`` (how the product leaves/accumulates) are forwarded to
    :func:`add_pathway_to_model`, so callers can dry-run on a copy and inspect the
    outcome without mutating ``host``."""
    return add_pathway_to_model(host.copy(), universal, reaction_ids,
                                match_namespace=match_namespace, rename=rename,
                                rename_map=rename_map, report=report, target=target,
                                ensure_exchange=ensure_exchange, gas_product=gas_product,
                                product_route=product_route)
