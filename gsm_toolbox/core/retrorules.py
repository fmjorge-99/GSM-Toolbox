"""Rule-based retrosynthesis (RetroRules-style) — the advanced engine that can propose
chemistry ABSENT from every database.

The classic engine (``pathway_search``) can only connect reactions that already exist
in a loaded database. That is why prodigiosin and paclitaxel fail even with BiGG +
MetaNetX + ModelSEED merged: the enzymes exist in nature but the *reactions* are not in
any of those databases. Reaction rules solve exactly this: a rule is a reaction pattern
(a reaction SMARTS at some atom "diameter"), and applying it to a target GENERATES a
hypothetical reaction and its precursors, rather than looking one up.

This module provides:

* :class:`ReactionRule` and :func:`load_rules` — read a RetroRules-format TSV, or the
  small curated ruleset bundled with the toolbox (:func:`bundled_rules`).
* :func:`apply_rule` — apply one rule retrosynthetically to a product, returning the
  precursor sets (RDKit does the SMARTS work).
* :func:`retro_expand` — a bounded retrosynthetic search: expand the target with the
  rules until every branch reaches a compound the host already has, or the step/'
  breadth budget is hit.
* :func:`assets_present` / :func:`full_ruleset_path` / :func:`download_hint` — the asset
  manager for the full RetroRules dataset (too large to bundle; fetched on demand).

Predictions are HYPOTHESES: a rule says "this transformation is chemically plausible",
not "this enzyme exists in your host". Callers must present rule-derived steps as
suggestions to verify, distinct from database-backed reactions.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

# RDKit is already a dependency (structure rendering). Import lazily so importing this
# module never fails on a stripped environment — the feature just reports unavailable.
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")          # rule SMARTS are noisy; we handle failures
    _RDKIT = True
except Exception:  # noqa: BLE001
    _RDKIT = False


def rdkit_available() -> bool:
    return _RDKIT


@dataclass(frozen=True)
class ReactionRule:
    """One retrosynthetic reaction rule.

    ``smarts`` is a reaction SMARTS in the RETRO direction — ``product>>precursors`` —
    so it is applied DIRECTLY to a target to yield its precursors (this matches how
    RetroRules stores its rules). ``diameter`` is the atom-environment radius (larger =
    more specific, fewer false hits). ``ec`` and ``name`` are informational.
    """

    rule_id: str
    smarts: str
    diameter: int = 0
    ec: str = ""
    name: str = ""
    score: float = 0.0     # RetroRules ``Score_normalized`` (0..1); higher = more reliable


# --- a small curated ruleset, bundled so the feature works with no download -----
# Each is a common biosynthetic transformation written in the RETRO direction
# (product>>precursors), so it is applied directly to a target to decompose it. These
# are deliberately generic (low diameter): enough to demonstrate rule-based expansion
# out of the box; the full RetroRules dataset (tens of thousands of rules) is fetched
# on demand for real coverage.
_BUNDLED: Tuple[Tuple[str, str, int, str, str], ...] = (
    ("BR_ester_hydrolysis", "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3].[O][C:4]",
     0, "3.1.-.-", "ester -> acid + alcohol"),
    ("BR_amide_hydrolysis", "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[O].[N:3]",
     0, "3.5.-.-", "amide -> acid + amine"),
    ("BR_carbonyl_reduction", "[C:1][O;H1]>>[C:1]=[O]", 0, "1.1.1.-",
     "alcohol <- aldehyde/ketone (reduction)"),
    ("BR_transamination", "[C:1]([NX3;H2])>>[C:1](=[O])", 0, "2.6.1.-",
     "amine <- keto (transamination)"),
    ("BR_o_demethylation", "[c,C:1][O:2][CH3]>>[c,C:1][O:2]", 0, "2.1.1.-",
     "methyl ether -> alcohol (demethylation)"),
    ("BR_dehydration", "[C:1]=[C:2]>>[C:1][C:2][O]", 0, "4.2.1.-",
     "alkene <- alcohol (dehydration)"),
)


def bundled_rules() -> List[ReactionRule]:
    """The curated ruleset shipped with the toolbox (works with no download)."""
    return [ReactionRule(rid, smarts, dia, ec, name)
            for rid, smarts, dia, ec, name in _BUNDLED]


def load_rules(path: str, *, max_diameter: Optional[int] = None,
               diameters: Optional[Sequence[int]] = None,
               limit: Optional[int] = None) -> List[ReactionRule]:
    """Load rules from a RetroRules-format TSV.

    RetroRules ships columns including ``# Rule_ID``, ``Rule_SMARTS``, ``Diameter``,
    ``Reaction_EC_number`` and ``Rule_usage``. We read defensively so minor schema
    differences between releases don't break loading, and skip rows without a usable
    SMARTS. ``diameters`` (a set) keeps only those atom-radii — the whole file is
    ~350k rules across 8 diameters, so filtering to one keeps the ruleset workable;
    ``limit`` caps the count.

    RetroRules stores each reaction in BOTH relative directions, so applying the rule's
    ``Rule_SMARTS`` DIRECTLY to a target (no reversal) covers retrosynthesis — matching
    how :func:`apply_rule` uses it and how the bundled rules are written.
    """
    want = set(int(d) for d in diameters) if diameters is not None else None
    rules: List[ReactionRule] = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower().lstrip("# ").strip(): c for c in (reader.fieldnames or [])}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_id = col("rule_id", "id")
        c_sm = col("rule_smarts", "reaction_smarts", "smarts")
        c_di = col("diameter")
        c_ec = col("reaction_ec_number", "ec_number", "ec", "ec_numbers")
        c_nm = col("rule_usage", "reaction_id", "name")
        c_sc = col("score_normalized", "score")
        if c_sm is None:
            raise ValueError("no SMARTS column found in the rules file")
        for row in reader:
            smarts = (row.get(c_sm) or "").strip()
            if not smarts or ">>" not in smarts:
                continue
            try:
                dia = int(row.get(c_di) or 0) if c_di else 0
            except ValueError:
                dia = 0
            if want is not None and dia not in want:
                continue
            if max_diameter is not None and dia > max_diameter:
                continue
            rid = (row.get(c_id) or f"rule{len(rules)}").strip() if c_id else \
                f"rule{len(rules)}"
            try:
                score = float(row.get(c_sc)) if c_sc and row.get(c_sc) else 0.0
            except ValueError:
                score = 0.0
            rules.append(ReactionRule(
                rid, smarts, dia,
                (row.get(c_ec) or "").strip() if c_ec else "",
                (row.get(c_nm) or "").strip() if c_nm else "",
                score))
            if limit is not None and len(rules) >= limit:
                break
    return rules


import re as _re

# Elements that appear on the reactant (left) side of a rule's SMARTS, cached per rule.
# Used to skip rules that cannot possibly match a target — a rule needing nitrogen is
# pointless on a target with no nitrogen. This prune is what makes searching the full
# 43k-rule diameter workable: it drops the great majority of rules for a given target
# before the expensive RDKit substructure match.
_RULE_ELEMENTS: Dict[str, frozenset] = {}
_ATOMIC = {1: "H", 6: "C", 7: "N", 8: "O", 15: "P", 16: "S", 9: "F",
           17: "Cl", 35: "Br", 53: "I", 5: "B", 14: "Si", 34: "Se"}


def _rule_reactant_elements(rule: "ReactionRule") -> frozenset:
    cached = _RULE_ELEMENTS.get(rule.rule_id)
    if cached is not None:
        return cached
    left = rule.smarts.split(">>", 1)[0]
    els = set()
    for num in _re.findall(r"#(\d+)", left):        # RetroRules writes [#6&v4...]
        els.add(int(num))
    if not els:                                     # organic-subset SMARTS (C, O, N…)
        for sym in _re.findall(r"[A-Z][a-z]?", left):
            els.add(sym)
    fs = frozenset(els)
    _RULE_ELEMENTS[rule.rule_id] = fs
    return fs


def _mol(smiles: str):
    if not _RDKIT:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:  # noqa: BLE001
        return None


def _mol_elements(mol) -> frozenset:
    """Atomic numbers present in a molecule (including implicit H)."""
    nums = set()
    hasH = False
    for a in mol.GetAtoms():
        nums.add(a.GetAtomicNum())
        if a.GetTotalNumHs() > 0:
            hasH = True
    if hasH:
        nums.add(1)
    return frozenset(nums)


def _rule_can_apply(rule: "ReactionRule", target_elements: frozenset) -> bool:
    """Cheap necessary condition: every element the rule's reactant needs must be in
    the target. Numeric (atomic-number) element sets compare directly; symbol-based
    fallbacks are not pruned (kept, to stay correct)."""
    need = _rule_reactant_elements(rule)
    if not need or not all(isinstance(e, int) for e in need):
        return True
    return need <= target_elements


def _canonical(smiles: str) -> str:
    m = _mol(smiles)
    return Chem.MolToSmiles(m) if m is not None else smiles


def _inchikey(smiles: str) -> str:
    m = _mol(smiles)
    if m is None:
        return ""
    try:
        return Chem.MolToInchiKey(m)
    except Exception:  # noqa: BLE001
        return ""


def apply_rule(rule: ReactionRule, product_smiles: str) -> List[List[str]]:
    """Apply one rule retrosynthetically to a product SMILES.

    Returns a list of precursor sets (each a list of SMILES). Empty if the rule does
    not match or RDKit is unavailable. The rule's SMARTS is already in the retro
    direction (product>>precursors), so it is applied directly.
    """
    if not _RDKIT:
        return []
    prod = _mol(product_smiles)
    if prod is None:
        return []
    try:
        rxn = AllChem.ReactionFromSmarts(rule.smarts)
        if rxn is None:
            return []
    except Exception:  # noqa: BLE001
        return []
    out: List[List[str]] = []
    seen: Set[Tuple[str, ...]] = set()
    try:
        for prodset in rxn.RunReactants((prod,)):
            smis = []
            ok = True
            for m in prodset:
                try:
                    Chem.SanitizeMol(m)
                    smis.append(Chem.MolToSmiles(m))
                except Exception:  # noqa: BLE001
                    ok = False
                    break
            if not ok or not smis:
                continue
            key = tuple(sorted(smis))
            if key in seen:
                continue
            seen.add(key)
            out.append(smis)
    except Exception:  # noqa: BLE001
        return out
    return out


@dataclass
class RuleRoute:
    """A rule-based retrosynthetic route to a target."""

    target: str = ""
    steps: List[dict] = field(default_factory=list)   # {rule_id, name, ec, product, precursors, score}
    terminal_precursors: List[str] = field(default_factory=list)  # native compounds it ends on
    complete: bool = False       # every branch reached a native compound

    @property
    def n_steps(self) -> int:
        return len(self.steps)

    @property
    def score(self) -> float:
        """Mean RetroRules reliability score over the route's steps (0 if unscored)."""
        scores = [s.get("score", 0.0) for s in self.steps]
        return sum(scores) / len(scores) if scores else 0.0

    @property
    def first_rule(self) -> str:
        """The rule that disconnects the target itself (defines the route's chemistry)."""
        for s in self.steps:
            if s.get("product") == self.target:
                return s.get("rule_id", "")
        return self.steps[-1].get("rule_id", "") if self.steps else ""


def ordered_rules(rules: Sequence[ReactionRule], seed: Optional[int] = None
                  ) -> List[ReactionRule]:
    """A DETERMINISTIC exploration order for a rule set (L1).

    The search can only sample a fraction of ~43k rules, so *which* rules it tries first
    decides which routes come back. Iterating a set/dict or relying on load order made the
    same query return different routes between runs, which made results impossible to
    reproduce or compare. Sorting by a stable key fixes the order; the ``seed`` then
    permutes it reproducibly, so a user can explore a different sample on purpose and
    still regenerate it exactly by re-using the seed.
    """
    # Primary key: higher reliability score first, then more specific (larger diameter),
    # then rule id — all stable across processes and platforms.
    base = sorted(rules, key=lambda r: (-getattr(r, "score", 0.0),
                                        -getattr(r, "diameter", 0),
                                        str(r.rule_id)))
    if not seed:
        return base
    import random
    rng = random.Random(int(seed))
    order = list(range(len(base)))
    rng.shuffle(order)
    return [base[i] for i in order]


def retro_expand(target_smiles: str, rules: Sequence[ReactionRule],
                 native_smiles: Optional[Sequence[str]] = None, *,
                 native_inchikeys: Optional[Sequence[str]] = None,
                 max_steps: int = 5, max_branch: int = 40,
                 time_budget: float = 120.0,
                 avoid_first: Optional[Set[str]] = None,
                 seed: Optional[int] = None) -> Optional[RuleRoute]:
    """Bounded rule-based retrosynthesis.

    Expand ``target_smiles`` with ``rules`` until every branch reaches a native compound
    or the depth/breadth budget is exhausted. Natives are given either as SMILES
    (``native_smiles``) or, for a real genome-scale model, directly as the host's
    ``inchi_key`` annotations (``native_inchikeys``) — cobra models usually carry
    InChIKeys but not SMILES. Returns the route (``complete=True`` if fully grounded in
    natives), or ``None`` if RDKit is unavailable or the target cannot be parsed.

    This is a foundation: a greedy depth-bounded search, not the full RetroRules
    scoring/beam search. It is enough to turn a target plus a ruleset into candidate
    heterologous steps, which is the capability the classic engine cannot provide.
    """
    if not _RDKIT:
        return None
    if _mol(target_smiles) is None:
        return None
    native_keys = set(native_inchikeys or [])
    native_keys |= {k for k in (_inchikey(s) for s in (native_smiles or [])) if k}
    # Match on the connectivity layer (first InChIKey block): stereochemistry and
    # protonation vary between a rule's output and a database entry for the same
    # compound, and would otherwise defeat exact-key matching.
    native_blocks = {k.split("-", 1)[0] for k in native_keys if k}

    import time
    deadline = time.monotonic() + max(1.0, time_budget)

    # Fixed exploration order — without this the same query returns different routes on
    # repeat runs (L1). Callers pass the user's configured seed.
    rules = ordered_rules(rules, seed)

    route = RuleRoute(target=_canonical(target_smiles))
    visited: Set[str] = set()

    def is_native(smiles: str) -> bool:
        k = _inchikey(smiles)
        return bool(k) and k.split("-", 1)[0] in native_blocks

    def expand(smiles: str, depth: int) -> bool:
        if is_native(smiles):
            if smiles not in route.terminal_precursors:
                route.terminal_precursors.append(smiles)
            return True
        if depth >= max_steps or time.monotonic() > deadline:
            return False           # depth or wall-clock budget hit — return what we have
        key = _canonical(smiles)
        if key in visited:
            return False
        visited.add(key)
        m = _mol(smiles)
        target_els = _mol_elements(m) if m is not None else frozenset()
        tried = 0
        for rule in rules:
            if time.monotonic() > deadline:
                return False
            # When enumerating alternative routes, skip the first-step rules already used
            # so each alternative disconnects the target with different chemistry.
            if depth == 0 and avoid_first and rule.rule_id in avoid_first:
                continue
            if not _rule_can_apply(rule, target_els):
                continue           # cheap prune: rule needs an element the target lacks
            # apply_rule's product ordering is not guaranteed stable across runs; sorting
            # each precursor set (and the sets themselves) keeps expansion reproducible.
            outcomes = sorted((tuple(sorted(p)) for p in apply_rule(rule, smiles)))
            for precursors in outcomes:
                tried += 1
                if tried > max_branch:
                    return False
                # Grounding all precursors reaches natives → accept this step.
                if all(is_native(p) or expand(p, depth + 1) for p in precursors):
                    route.steps.append({
                        "rule_id": rule.rule_id, "name": rule.name, "ec": rule.ec,
                        "product": key, "precursors": [_canonical(p) for p in precursors],
                        "score": getattr(rule, "score", 0.0),
                    })
                    for p in precursors:
                        if is_native(p) and p not in route.terminal_precursors:
                            route.terminal_precursors.append(p)
                    return True
        return False

    route.complete = expand(target_smiles, 0)
    route.steps.reverse()          # precursor-first, like the classic engine
    return route


# --- plausibility filtering (VI.3A) ---------------------------------------------------
# Precursors that cannot be the carbon source of an organic product. A rule route that
# "grounds" a carbon skeleton in nitrite or O2 has not found chemistry — it has found an
# artefact of applying a rule outside its domain.
_INORGANIC_PRECURSORS = {
    "O=N[O-]", "[O-][N+](=O)[O-]", "O=[N+]([O-])[O-]", "N", "O=O", "O", "[H]O[H]",
    "O=C=O", "[C-]#[O+]", "S", "P", "[NH4+]", "N#N",
}


def _n_carbons(smiles: str) -> int:
    """Carbon count of a SMILES (0 if unparseable)."""
    try:
        m = _mol(smiles)
        if m is None:
            return 0
        return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
    except Exception:  # noqa: BLE001
        return 0


def _is_inorganic(smiles: str) -> bool:
    if not smiles:
        return True
    if smiles in _INORGANIC_PRECURSORS:
        return True
    return _n_carbons(smiles) == 0


def route_plausibility(route: "RuleRoute") -> Tuple[bool, str]:
    """Is this rule route chemically credible? Returns ``(plausible, reason)``.

    Three cheap, decisive tests, drawn from the failure modes actually observed:

    1. **Inorganic terminal precursor** — a carbon skeleton cannot come from nitrite,
       nitrate or O2. The default-seed n-butanol route grounded in nitrite.
    2. **Carbon must not be created from nothing** — a step whose product has more
       carbons than all its precursors combined is impossible.
    3. **Gross carbon mismatch** — a step that discards most of the carbon it was given
       (e.g. C18 → C4 in one "synthesis" step) is almost always a misapplied rule.
    """
    if route is None or not route.steps:
        return False, "no steps"
    for p in route.terminal_precursors:
        if _is_inorganic(p):
            return False, (f"grounds the carbon skeleton in an inorganic precursor "
                           f"({p}) — a rule applied outside its chemical domain")
    for s in route.steps:
        prod_c = _n_carbons(s.get("product", ""))
        prec_c = [_n_carbons(p) for p in s.get("precursors", [])]
        total_prec = sum(prec_c)
        if prod_c <= 0:
            continue
        if total_prec == 0:
            return False, (f"a step makes a C{prod_c} product from precursors with no "
                           "carbon at all")
        if prod_c > total_prec:
            return False, (f"a step creates carbon from nothing (C{total_prec} → "
                           f"C{prod_c})")
        biggest = max(prec_c) if prec_c else 0
        if biggest >= 2 * max(1, prod_c) and biggest - prod_c >= 4:
            return False, (f"a step discards most of its carbon (C{biggest} → "
                           f"C{prod_c}), which is not a synthesis step")
    return True, ""


def filter_plausible(routes: Sequence["RuleRoute"]) -> Tuple[List["RuleRoute"], List[str]]:
    """Split routes into the credible ones and the reasons the rest were rejected."""
    keep, rejected = [], []
    for r in routes:
        ok, why = route_plausibility(r)
        if ok:
            keep.append(r)
        else:
            rejected.append(why)
    return keep, rejected


# Ranking keys for alternative routes; each returns a sort tuple (smaller = better).
RANK_KEYS = {
    "score": lambda r: (-r.score, r.n_steps),                 # most reliable rules first
    "steps": lambda r: (r.n_steps, -r.score),                 # shortest route first
    "precursors": lambda r: (len(r.terminal_precursors), r.n_steps, -r.score),
}
RANK_LABELS = {
    "score": "RetroRules reliability score",
    "steps": "fewest steps",
    "precursors": "fewest native precursors",
}


def rank_routes(routes: Sequence["RuleRoute"], by: str = "score") -> List["RuleRoute"]:
    key = RANK_KEYS.get(by, RANK_KEYS["score"])
    return sorted(routes, key=key)


def retro_expand_multi(target_smiles: str, rules: Sequence[ReactionRule],
                       *, n_alternatives: int = 3, rank_by: str = "score",
                       native_smiles: Optional[Sequence[str]] = None,
                       native_inchikeys: Optional[Sequence[str]] = None,
                       max_steps: int = 5, max_branch: int = 40,
                       time_budget: float = 120.0,
                       seed: Optional[int] = None) -> List[RuleRoute]:
    """Enumerate up to ``n_alternatives`` DISTINCT rule-based routes and rank them.

    Each alternative is forced to disconnect the target with a different first-step rule,
    so the results are genuinely different chemistries rather than trivial variants. The
    per-route time budget is shared across attempts. Ranking is by ``rank_by`` (see
    ``RANK_KEYS``): reliability score, fewest steps, or fewest native precursors.
    """
    import time
    if not _RDKIT or _mol(target_smiles) is None:
        return []
    deadline = time.monotonic() + max(1.0, time_budget)
    found: List[RuleRoute] = []
    used_first: Set[str] = set()
    seen_signatures: Set[tuple] = set()
    # Try a few more times than requested: some attempts yield an already-seen or empty
    # route, and we want n genuinely distinct ones if they exist.
    for _ in range(max(1, n_alternatives) + 4):
        if time.monotonic() > deadline or len(found) >= max(1, n_alternatives):
            break
        remaining = max(1.0, deadline - time.monotonic())
        route = retro_expand(target_smiles, rules, native_smiles=native_smiles,
                             native_inchikeys=native_inchikeys, max_steps=max_steps,
                             max_branch=max_branch, time_budget=remaining,
                             avoid_first=used_first, seed=seed)
        if route is None or not route.steps:
            break
        sig = tuple(sorted(s.get("rule_id", "") for s in route.steps))
        if sig not in seen_signatures:
            seen_signatures.add(sig)
            found.append(route)
        first = route.first_rule
        if first:
            used_first.add(first)
        else:
            break                  # can't diversify further without a first-rule handle
    return rank_routes(found, rank_by)


# --- asset manager for the full RetroRules dataset ------------------------------
# The full RetroRules RR02 ruleset is ~350k rules (502 MB uncompressed), too large to
# bundle, so it is downloaded on demand into the user's data dir. A missing dataset
# disables ONLY the full-coverage rules — the bundled curated set still works.

# RetroRules RR02, non-stereo, diameters 2-16 (Zenodo, ~43 MB compressed).
RETRORULES_URL = ("https://zenodo.org/records/5827969/files/"
                  "retrorules_rr02_rp3_nohs.tar.gz?download=1")
RETRORULES_MD5 = "a881165602bb9e3c416013df3cb3dce1"
DEFAULT_DIAMETER = 6      # a workable balance of specificity vs generality

# Bump whenever the FIELDS parsed into ReactionRule change, so stale per-diameter pickles
# are rebuilt rather than silently supplying rules missing the new data (L12).
# v2: added ReactionRule.score (RetroRules Score_normalized).
RULE_CACHE_VERSION = 2


def rules_dir() -> str:
    from .cache import databases_dir
    d = os.path.join(os.path.dirname(databases_dir()), "retrorules")
    os.makedirs(d, exist_ok=True)
    return d


def full_ruleset_path() -> str:
    """The extracted flat TSV of all RetroRules rules (present once downloaded)."""
    return os.path.join(rules_dir(), "retrorules_rr02_rp3_nohs",
                        "retrorules_rr02_flat_all.tsv")


def archive_path() -> str:
    return os.path.join(rules_dir(), "retrorules_rr02_rp3_nohs.tar.gz")


def assets_present() -> bool:
    p = full_ruleset_path()
    return os.path.exists(p) and os.path.getsize(p) > 0


def download_hint() -> str:
    return ("The full RetroRules dataset is not installed. It can be downloaded "
            f"(~43 MB) from Zenodo:\n  {RETRORULES_URL}\n"
            f"and extracted into:\n  {rules_dir()}\n"
            "Or call retrorules.install_full_ruleset(). The bundled curated ruleset "
            "works without it, with limited coverage.")


def install_full_ruleset(progress=None) -> str:
    """Download + verify + extract the full RetroRules dataset. Returns the TSV path.

    ``progress`` is an optional callable(message:str) for UI feedback. Idempotent: if
    the dataset is already extracted, returns immediately.
    """
    import tarfile
    import urllib.request

    if assets_present():
        return full_ruleset_path()

    def say(m):
        if progress:
            try:
                progress(m)
            except Exception:  # noqa: BLE001
                pass

    arc = archive_path()
    if not (os.path.exists(arc) and _md5(arc) == RETRORULES_MD5):
        say("Downloading RetroRules dataset (~43 MB)…")
        urllib.request.urlretrieve(RETRORULES_URL, arc)
        digest = _md5(arc)
        if digest != RETRORULES_MD5:
            raise ValueError(
                f"RetroRules download is corrupt (md5 {digest}, expected "
                f"{RETRORULES_MD5}). Delete {arc} and try again.")
    say("Extracting rules…")
    with tarfile.open(arc, "r:gz") as tf:
        tf.extractall(rules_dir())
    if not assets_present():
        raise ValueError("RetroRules archive extracted but the rules TSV is missing.")
    return full_ruleset_path()


def _md5(path: str) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_full_ruleset(*, diameter: int = DEFAULT_DIAMETER,
                      limit: Optional[int] = None) -> List[ReactionRule]:
    """Load the full dataset filtered to a single ``diameter`` (~43k rules per
    diameter), caching a compact per-diameter pickle so later loads are fast (parsing
    the 502 MB TSV every time would be wasteful)."""
    import pickle
    # The cache filename carries RULE_CACHE_VERSION: when the parsed fields change (as
    # they did when Score_normalized was added), an old pickle would silently supply
    # rules with score 0 and the "rank by reliability score" option would degrade to
    # step-count with no warning. Bumping the version forces a rebuild instead. (L12)
    cache = os.path.join(rules_dir(), f"rr02_d{diameter}_v{RULE_CACHE_VERSION}.pkl")
    if limit is None and os.path.exists(cache):
        try:
            with open(cache, "rb") as fh:
                payload = pickle.load(fh)
            if isinstance(payload, dict) and payload.get("version") == RULE_CACHE_VERSION:
                return payload["rules"]
        except Exception:  # noqa: BLE001
            pass
    rules = load_rules(full_ruleset_path(), diameters=[diameter], limit=limit)
    if limit is None:
        try:
            with open(cache, "wb") as fh:
                pickle.dump({"version": RULE_CACHE_VERSION, "diameter": diameter,
                             "rules": rules}, fh)
        except Exception:  # noqa: BLE001
            pass
        _purge_stale_rule_caches(diameter)
    return rules


def _purge_stale_rule_caches(diameter: int) -> None:
    """Delete rule caches written by an older schema so they can't be picked up again
    and don't waste hundreds of MB on disk."""
    try:
        d = rules_dir()
        keep = f"rr02_d{diameter}_v{RULE_CACHE_VERSION}.pkl"
        for f in os.listdir(d):
            if f.startswith(f"rr02_d{diameter}") and f.endswith(".pkl") and f != keep:
                try:
                    os.remove(os.path.join(d, f))
                except OSError:
                    pass
    except Exception:  # noqa: BLE001
        pass


def active_rules(*, diameter: int = DEFAULT_DIAMETER,
                 limit: Optional[int] = None) -> List[ReactionRule]:
    """The rules to search with: the full dataset (one diameter) if installed, else the
    bundled curated set."""
    if assets_present():
        try:
            return load_full_ruleset(diameter=diameter, limit=limit)
        except Exception:  # noqa: BLE001
            pass
    return bundled_rules()


# --- host-facing integration ----------------------------------------------------

def host_inchikeys(model) -> List[str]:
    """Every InChIKey the host model carries (cobra models annotate ``inchi_key``).

    These are the compounds a rule-based route may bottom out in — the same "native
    precursors" idea as the classic engine, matched structurally rather than by id.
    """
    keys: List[str] = []
    for m in model.metabolites:
        ann = getattr(m, "annotation", None) or {}
        v = ann.get("inchi_key") or ann.get("inchikey") or ann.get("InChIKey")
        if not v:
            continue
        for k in (v if isinstance(v, (list, tuple)) else [v]):
            k = str(k).strip()
            if k:
                keys.append(k)
    return keys


def metabolite_smiles(met) -> str:
    """A SMILES for a cobra metabolite, from its SMILES or InChI annotation (many
    models carry one or the other). Empty if none is available or RDKit is absent."""
    ann = getattr(met, "annotation", None) or {}

    def first(*names):
        for n in names:
            v = ann.get(n)
            if v:
                return (v[0] if isinstance(v, (list, tuple)) else str(v)).strip()
        return ""

    smi = first("smiles", "SMILES")
    if smi:
        return smi
    if _RDKIT:
        inchi = first("inchi", "InChI")
        if inchi:
            try:
                m = Chem.MolFromInchi(inchi)
                if m is not None:
                    return Chem.MolToSmiles(m)
            except Exception:  # noqa: BLE001
                pass
    return ""


def _formula(smiles: str) -> str:
    if not _RDKIT:
        return ""
    m = _mol(smiles)
    if m is None:
        return ""
    try:
        from rdkit.Chem import rdMolDescriptors
        return rdMolDescriptors.CalcMolFormula(m)
    except Exception:  # noqa: BLE001
        return ""


def _smiles_slug(smiles: str, seen: Dict[str, str]) -> str:
    """A short, stable, unique metabolite id for a SMILES (InChIKey-based)."""
    if smiles in seen:
        return seen[smiles]
    key = _inchikey(smiles)
    base = ("rr_" + key.split("-", 1)[0]).lower() if key else f"rr_cpd{len(seen)}"
    mid = f"{base}_c"
    i = 2
    while mid in seen.values():
        mid = f"{base}_{i}_c"
        i += 1
    seen[smiles] = mid
    return mid


# Words that carry no information about the chemistry, so they must never become an id
# on their own — "RR_retro" told a reader nothing (VI.4).
_UNINFORMATIVE_SLUGS = {"retro", "retrorules", "step", "rule", "reaction", "synthesis",
                        "unknown", "product"}


def _readable_rxn_id(name: str, index: int, used: "Set[str]", *,
                     rule_id: str = "", product_smiles: str = "") -> str:
    """A traceable reaction id derived from what the step actually does (L13, VI.4).

    Preference order, so the id always says something:

    1. the reaction's human name — ``"Tryptamine synthesis"`` → ``RR_tryptamine_synthesis``
    2. the product's molecular formula — ``RR_C4H10O_synthesis``
    3. the RetroRules rule id — ``RR_rule_02_915c8ea59c``
    4. only then the step counter.
    """
    slug = _re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip()).strip("_").lower()[:44]
    if slug and slug not in _UNINFORMATIVE_SLUGS:
        base = f"RR_{slug}"
    else:
        formula = _formula(product_smiles) if product_smiles else ""
        if formula:
            base = f"RR_{_re.sub(r'[^A-Za-z0-9]+', '', formula)}_synthesis"
        elif rule_id:
            base = "RR_rule_" + _re.sub(r"[^A-Za-z0-9]+", "_",
                                        str(rule_id)).strip("_").lower()[:40]
        else:
            base = f"RR_step{index}"
    rid = base
    n = 2
    while rid in used:
        rid = f"{base}_{n}"
        n += 1
    return rid


def build_suggested_model(route: "RuleRoute", target_smiles: str, *,
                          target_name: str = "", names=None):
    """Turn a rule-based :class:`RuleRoute` into a cobra model + reaction ids, so the
    normal *Apply pathway* flow can add it to the host.

    Each unique SMILES becomes a metabolite (id from its InChIKey, formula from RDKit,
    an ``inchi_key`` annotation so `apply_pathway` can reconcile it with host natives);
    each step becomes a reaction (precursors → product). ``names`` optionally maps a
    SMILES to a human-readable name. Returns ``(model, reaction_ids, target_met_id)``.
    """
    import cobra

    names = names or {}
    model = cobra.Model("retrorules_suggested")
    seen: Dict[str, str] = {}
    mets: Dict[str, cobra.Metabolite] = {}

    def met_for(smiles: str) -> cobra.Metabolite:
        mid = _smiles_slug(smiles, seen)
        if mid not in mets:
            key = _inchikey(smiles)
            m = cobra.Metabolite(mid, name=names.get(smiles) or _canonical(smiles),
                                 formula=_formula(smiles), compartment="c")
            if key:
                m.annotation["inchi_key"] = key
            m.annotation["smiles"] = _canonical(smiles)
            mets[mid] = m
        return mets[mid]

    reaction_ids: List[str] = []
    used_ids: Set[str] = set()
    for i, s in enumerate(route.steps, 1):
        prod = met_for(s["product"])
        precs = [met_for(p) for p in s["precursors"]]
        # A readable reaction name: prefer a name built from the product, so it is
        # meaningful (and becomes the default suggested id). Fall back to EC / rule id.
        prod_name = names.get(s["product"]) or (
            target_name if s["product"] == target_smiles else "")
        if prod_name:
            rname = f"{prod_name} synthesis"
        elif s.get("name"):
            rname = s["name"]
        elif s.get("ec"):
            rname = f"EC {s['ec']}"
        else:
            rname = f"RetroRules step {i}"
        # The ID follows the NAME, not a bare counter: "RR_step2" tells a reader nothing,
        # whereas "RR_tryptamine_synthesis" is traceable in a model, a table or a report
        # (L13). Uniqueness is enforced with a numeric suffix only when needed.
        rid = _readable_rxn_id(rname, i, used_ids, rule_id=s.get("rule_id", ""),
                               product_smiles=s.get("product", ""))
        used_ids.add(rid)
        rxn = cobra.Reaction(rid, name=rname)
        rxn.lower_bound, rxn.upper_bound = 0.0, 1000.0
        coeffs = {p: -1.0 for p in precs}
        coeffs[prod] = coeffs.get(prod, 0.0) + 1.0
        rxn.add_metabolites(coeffs)
        rxn.annotation["retrorules.rule"] = s.get("rule_id", "")
        if s.get("ec"):
            rxn.annotation["ec-code"] = s["ec"]
        model.add_reactions([rxn])
        reaction_ids.append(rid)

    target_mid = _smiles_slug(target_smiles, seen) if target_smiles in seen else ""
    if not target_mid and target_smiles:
        # target may equal a product SMILES already created
        tkey = _inchikey(target_smiles)
        for mid, m in mets.items():
            if m.annotation.get("inchi_key") == tkey:
                target_mid = mid
                break
    if target_name and target_mid and target_mid in mets:
        mets[target_mid].name = target_name
    return model, reaction_ids, target_mid


def suggest_routes(target_smiles: str, host, *, diameter: int = DEFAULT_DIAMETER,
                   max_steps: int = 4, max_branch: int = 60, time_budget: float = 120.0,
                   rules: Optional[Sequence[ReactionRule]] = None) -> Optional[RuleRoute]:
    """Rule-based retrosynthetic suggestion for a target, grounded in a host model.

    Ties the pieces together: pull the host's native InChIKeys, load the active
    ruleset, and expand. Returns a :class:`RuleRoute` (``complete`` if every branch
    reached a native compound), or ``None`` if RDKit is unavailable or the target
    SMILES cannot be parsed. Rule-derived steps are HYPOTHESES to verify, not
    database-backed reactions.
    """
    if not _RDKIT:
        return None
    rules = rules if rules is not None else active_rules(diameter=diameter)
    return retro_expand(target_smiles, rules,
                        native_inchikeys=host_inchikeys(host),
                        max_steps=max_steps, max_branch=max_branch,
                        time_budget=time_budget)


def suggest_routes_multi(target_smiles: str, host, *, diameter: int = DEFAULT_DIAMETER,
                         n_alternatives: int = 3, rank_by: str = "score",
                         max_steps: int = 4, max_branch: int = 60,
                         time_budget: float = 120.0,
                         rules: Optional[Sequence[ReactionRule]] = None,
                         seed: Optional[int] = None
                         ) -> List[RuleRoute]:
    """Several ranked rule-based routes to a target, grounded in a host model.

    Like :func:`suggest_routes` but returns up to ``n_alternatives`` genuinely different
    routes, ranked by ``rank_by`` (``"score"`` / ``"steps"`` / ``"precursors"``). Empty
    list if RDKit is unavailable or the target SMILES cannot be parsed.
    """
    if not _RDKIT:
        return []
    if seed is None:                    # default to the user's configured seed (L1)
        try:
            from . import preferences
            seed = preferences.retrorules_seed()
        except Exception:  # noqa: BLE001
            seed = None
    rules = rules if rules is not None else active_rules(diameter=diameter)
    routes = retro_expand_multi(target_smiles, rules, n_alternatives=n_alternatives,
                                rank_by=rank_by, native_inchikeys=host_inchikeys(host),
                                max_steps=max_steps, max_branch=max_branch,
                                time_budget=time_budget, seed=seed)
    keep, _rejected = filter_plausible(routes)
    return keep


def suggest_routes_multiseed(target_smiles: str, host, *, seeds: Sequence[int],
                             diameter: int = DEFAULT_DIAMETER,
                             n_alternatives: int = 3, rank_by: str = "score",
                             max_steps: int = 4, max_branch: int = 60,
                             time_budget_per_seed: float = 90.0,
                             rules: Optional[Sequence[ReactionRule]] = None,
                             progress: Optional[callable] = None) -> dict:
    """Search with several seeds and return only the chemically plausible routes (VI.3A).

    The rule search can sample only a fraction of ~43k rules, so which routes come back
    depends on the exploration order — one seed returned a route grounding n-butanol in
    nitrite while another found the textbook butyrate route at a *higher* score. Running
    a handful of seeds and pooling the survivors is therefore the honest default.

    Returns ``{"routes": [...], "rejected": [reasons], "seeds": [...], "by_seed": {...}}``
    with routes de-duplicated by their rule signature and ranked by ``rank_by``.
    """
    if not _RDKIT:
        return {"routes": [], "rejected": [], "seeds": list(seeds), "by_seed": {}}
    rules = rules if rules is not None else active_rules(diameter=diameter)
    natives = host_inchikeys(host)

    pooled: List[RuleRoute] = []
    seen_sig: Set[tuple] = set()
    rejected: List[str] = []
    by_seed: Dict[int, int] = {}
    for i, sd in enumerate(seeds, 1):
        if progress:
            progress(f"Searching reaction rules — seed {sd} ({i}/{len(seeds)})…")
        found = retro_expand_multi(target_smiles, rules, n_alternatives=n_alternatives,
                                   rank_by=rank_by, native_inchikeys=natives,
                                   max_steps=max_steps, max_branch=max_branch,
                                   time_budget=time_budget_per_seed, seed=sd)
        keep, why = filter_plausible(found)
        rejected.extend(why)
        by_seed[sd] = len(keep)
        for r in keep:
            sig = tuple(sorted(s.get("rule_id", "") for s in r.steps))
            if sig in seen_sig:
                continue
            seen_sig.add(sig)
            r.seed = sd                 # remember which seed produced it
            pooled.append(r)
    return {"routes": rank_routes(pooled, rank_by), "rejected": rejected,
            "seeds": list(seeds), "by_seed": by_seed}
