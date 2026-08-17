"""Propose correspondences between host metabolites and database metabolites.

The toolbox already unifies compounds across namespaces in two places: the union-find in
:class:`pathway_search.Canonicalizer`, and :func:`namespace.match_database_to_host`. Both
start from **cross-reference tokens** — a candidate pair only exists if the two entries
already share a KEGG, ChEBI, MetaNetX, SEED or BiGG identifier. That covers the great
majority of cases and this module deliberately does not duplicate it.

What it adds is the case those two cannot reach: **no shared cross-reference at all**.

The worked example is ferredoxin, and it is worth stating because it shaped every design
decision here. Routing fucoxanthin through *Synechocystis* iJN678 needs one database
reaction that consumes reduced ferredoxin. The host calls it ``fdxrd_c`` and annotates it
with BioCyc and ChEBI; the merged database calls it ``Reduced_ferredoxins_c`` and
annotates it with MetaNetX and SEED. The two annotation sets are **disjoint**, so no
candidate pair is ever generated, the metabolite is added to the model as an orphan, and
the whole route carries zero flux while looking perfectly well-formed.

Worse, the two obvious fallbacks both fail on it:

* **Formula equality** says *no*: the host writes ``Fe2S2X`` and the database writes
  ``Fe2R8S2``. Neither is a real molecular formula — ferredoxin is a protein-bound
  cofactor and both are placeholders with an R/X group standing for the apoprotein. A
  rule of "never merge across a formula difference" rejects precisely the match that
  matters.
* **Name containment** is unsafe in the other direction: ``Zeaxanthin`` is a substring of
  ``Zeaxanthin diglucoside``, which is a different compound.

So the design is: several independent kinds of evidence, each scored; full normalised
name equality rather than containment; formula treated as *supporting* evidence that can
veto, except for a curated list of protein-bound cofactors where the formula is known to
be a placeholder; and nothing applied automatically below a high confidence threshold.

Nothing here mutates either model. A proposal is a suggestion about two existing objects;
the accepted mapping is a separate artefact the caller stores and can revoke.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import cobra

from . import namespace

# --------------------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------------------
#: A pair at or above this is safe to apply without asking. It takes either a full
#: InChIKey match, or a cross-reference match backed by agreeing chemistry — never a name
#: on its own, however exact.
AUTO_ACCEPT = 90.0

#: Below this a pair is not worth showing; the list would fill with noise.
REVIEW_FLOOR = 35.0

_W_INCHIKEY_FULL = 100.0     # same stereochemistry and skeleton: decisive
_W_INCHIKEY_BLOCK = 60.0     # same skeleton, possibly different stereoisomer
_W_XREF_EACH = 25.0          # per shared cross-reference
_W_XREF_CAP = 75.0           # three agreeing references is as convincing as it gets
_W_FORMULA = 30.0            # supporting only: C40H56 fits many carotenoids
_W_NAME = 40.0               # full normalised equality, never containment
_W_COFACTOR = 35.0           # both sides are the same named protein-bound cofactor

_P_INCHIKEY_CONFLICT = -80.0  # different skeletons — almost certainly different compounds
_P_FORMULA_CONFLICT = -60.0   # real formulas that disagree

#: Protein-bound cofactors whose "formula" is a placeholder standing for the apoprotein,
#: so elemental comparison between two databases is meaningless. Matching one of these
#: names on both sides suspends the formula veto; it does not by itself accept the pair.
#:
#: Curated deliberately rather than inferred. A heuristic such as "contains an R group"
#: would also catch genuine R-group chemistry — acyl-CoA templates, lipid classes — where
#: a formula difference really does mean different compounds.
COFACTOR_CLASSES: Dict[str, Tuple[str, ...]] = {
    "ferredoxin": ("ferredoxin", "ferredoxins"),
    "flavodoxin": ("flavodoxin", "flavodoxins"),
    "thioredoxin": ("thioredoxin", "thioredoxins"),
    "glutaredoxin": ("glutaredoxin", "glutaredoxins"),
    "acyl-carrier-protein": ("acyl carrier protein", "acyl-carrier protein", "acp"),
    "plastocyanin": ("plastocyanin", "plastocyanins"),
    "cytochrome-c": ("cytochrome c", "cytochrome-c", "ferrocytochrome c",
                     "ferricytochrome c"),
    "rubredoxin": ("rubredoxin", "rubredoxins"),
}

#: Redox state must still agree — reduced ferredoxin is not oxidised ferredoxin. These
#: are checked as whole words so "oxidised" inside a longer name cannot be missed.
_REDOX_TERMS = (("reduced", "red"), ("oxidized", "oxidised", "ox", "oxid"))

#: Greek letters spelled out, as a name may use either form or the symbol.
_GREEK = {
    "alpha": "a", "α": "a", "beta": "b", "β": "b", "gamma": "g", "γ": "g",
    "delta": "d", "δ": "d", "epsilon": "e", "ε": "e", "zeta": "z", "ζ": "z",
    "omega": "w", "ω": "w", "psi": "psi", "ψ": "psi",
}

_COMPARTMENT_SUFFIX = re.compile(r"_[a-z][a-z0-9]?$", re.I)


# --------------------------------------------------------------------------------------
# Evidence and proposals
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Evidence:
    """One reason to believe (or disbelieve) a pair, and what it contributed."""

    kind: str            # xref | inchikey | inchikey_block | formula | name | cofactor
    detail: str          # human-readable, shown verbatim in the review dialog
    weight: float

    def __str__(self) -> str:
        return f"{self.detail} ({self.weight:+.0f})"


@dataclass
class Proposal:
    """A suggested correspondence between one host and one database metabolite."""

    host_id: str
    db_id: str
    score: float
    evidence: List[Evidence] = field(default_factory=list)
    #: Set when something disqualifies the pair from automatic acceptance. The proposal
    #: is still shown — a blocked pair the user recognises as correct is exactly the case
    #: an override exists for — but it can never be applied without one.
    blocked: str = ""
    host_name: str = ""
    db_name: str = ""
    compartment: str = ""

    @property
    def confidence(self) -> str:
        if self.blocked:
            return "blocked"
        if self.score >= AUTO_ACCEPT:
            return "high"
        if self.score >= (AUTO_ACCEPT + REVIEW_FLOOR) / 2:
            return "medium"
        return "low"

    @property
    def auto_acceptable(self) -> bool:
        return not self.blocked and self.score >= AUTO_ACCEPT

    def why(self) -> str:
        """One line naming every piece of evidence, for the dialog and for reports."""
        parts = [str(e) for e in sorted(self.evidence, key=lambda e: -abs(e.weight))]
        text = "; ".join(parts) or "no evidence"
        return f"{text} — BLOCKED: {self.blocked}" if self.blocked else text


# --------------------------------------------------------------------------------------
# Normalisation
# --------------------------------------------------------------------------------------
def normalise_name(text: str) -> str:
    """Casefold a compound name to a comparable form.

    Strips a trailing compartment suffix, expands spelled-out Greek prefixes so
    ``beta_Carotene``, ``b-carotene`` and ``β-carotene`` agree, and collapses every
    separator. Deliberately *not* stemming or removing words: dropping "diglucoside"
    would make ``Zeaxanthin diglucoside`` equal to ``Zeaxanthin``.
    """
    if not text:
        return ""
    out = _COMPARTMENT_SUFFIX.sub("", str(text).strip())
    out = out.casefold()
    for word, letter in _GREEK.items():
        # Only as a standalone prefix token, so "betaine" is not rewritten to "baine".
        out = re.sub(rf"(?:^|(?<=[\s_\-,\(]))\{word}\b" if word.startswith("\\")
                     else rf"(?:^|(?<=[\s_\-,\(])){re.escape(word)}(?=[\s_\-]|$)",
                     letter, out)
    out = re.sub(r"[\s_\-,'\.\(\)\[\]]+", "", out)
    return out


def _plural_insensitive(name: str) -> str:
    """``ferredoxins`` and ``ferredoxin`` are the same cofactor pool."""
    return name[:-1] if name.endswith("s") and len(name) > 3 else name


def _formula_key(met) -> Optional[tuple]:
    return namespace._formula_key(met)


def _is_placeholder_formula(met) -> bool:
    """True when the formula stands for a macromolecule rather than a molecule.

    Both ``Fe2S2X`` and ``Fe2R8S2`` describe ferredoxin; the R/X is the apoprotein. Two
    databases will not agree on such a string, so comparing them elementally is noise.
    """
    formula = (getattr(met, "formula", "") or "")
    return bool(re.search(r"[RX](?![a-z])", formula))


def _inchikey_full(met) -> Optional[str]:
    ann = getattr(met, "annotation", None) or {}
    for key in ("inchi_key", "inchikey", "InChI Key", "inchikey_full"):
        value = ann.get(key)
        if isinstance(value, (list, tuple)):
            value = value[0] if value else None
        if value:
            text = str(value).strip().split("/")[-1].upper()
            return text if text.count("-") >= 1 else None
    return None


def cofactor_class(met) -> Optional[str]:
    """Which curated cofactor class this metabolite belongs to, if any."""
    haystack = f"{getattr(met, 'name', '') or ''} {getattr(met, 'id', '') or ''}".casefold()
    haystack = re.sub(r"[\s_\-]+", " ", haystack)
    for label, terms in COFACTOR_CLASSES.items():
        if any(re.search(rf"\b{re.escape(t)}\b", haystack) for t in terms):
            return label
    return None


def _redox_state(met) -> Optional[str]:
    """'reduced' / 'oxidized' / None, so the two halves of a couple never merge."""
    text = re.sub(r"[\s_\-]+", " ",
                  f"{getattr(met, 'name', '') or ''} {getattr(met, 'id', '') or ''}"
                  .casefold())
    for state, terms in zip(("reduced", "oxidized"), _REDOX_TERMS):
        if any(re.search(rf"\b{t}\b", text) for t in terms):
            return state
    return None


def _compartments_agree(a, b) -> bool:
    ca = (getattr(a, "compartment", "") or "").casefold()
    cb = (getattr(b, "compartment", "") or "").casefold()
    if not ca or not cb:
        return True                     # an unstated compartment is not a disagreement
    return ca == cb


# --------------------------------------------------------------------------------------
# Scoring one pair
# --------------------------------------------------------------------------------------
def score_pair(host_met, db_met) -> Proposal:
    """Weigh every kind of evidence for one candidate pair.

    Public, because it is the unit the tests pin down and the dialog explains.
    """
    proposal = Proposal(
        host_id=host_met.id, db_id=db_met.id, score=0.0,
        host_name=getattr(host_met, "name", "") or "",
        db_name=getattr(db_met, "name", "") or "",
        compartment=(getattr(host_met, "compartment", "") or ""))

    if not _compartments_agree(host_met, db_met):
        proposal.blocked = "different compartments"
        return proposal

    add = proposal.evidence.append

    # --- cross-references -------------------------------------------------------------
    shared = (set(namespace.metabolite_tokens(host_met))
              & set(namespace.metabolite_tokens(db_met)))
    # A shared token derived from the id itself ("bigg:<id>") is not independent evidence
    # when the ids are literally equal; it is the thing being questioned.
    informative = {t for t in shared if not t.startswith("base:")}
    if informative:
        weight = min(_W_XREF_EACH * len(informative), _W_XREF_CAP)
        add(Evidence("xref", f"{len(informative)} shared reference(s): "
                             f"{', '.join(sorted(informative)[:3])}", weight))
        proposal.score += weight

    # --- structure --------------------------------------------------------------------
    hk, dk = _inchikey_full(host_met), _inchikey_full(db_met)
    hb, db_block = namespace._inchikey_block(host_met), namespace._inchikey_block(db_met)
    if hk and dk and hk == dk:
        add(Evidence("inchikey", f"identical InChIKey {hk}", _W_INCHIKEY_FULL))
        proposal.score += _W_INCHIKEY_FULL
    elif hb and db_block:
        if hb == db_block:
            add(Evidence("inchikey_block",
                         f"same skeleton {hb} (stereochemistry may differ)",
                         _W_INCHIKEY_BLOCK))
            proposal.score += _W_INCHIKEY_BLOCK
        else:
            # This is the guard that keeps two different C40H56 carotenoids apart when
            # their formulas — and even some of their references — agree.
            add(Evidence("inchikey_block",
                         f"different skeletons {hb} vs {db_block}",
                         _P_INCHIKEY_CONFLICT))
            proposal.score += _P_INCHIKEY_CONFLICT
            proposal.blocked = "different chemical skeletons (InChIKey)"

    # --- formula and charge -----------------------------------------------------------
    placeholder = _is_placeholder_formula(host_met) or _is_placeholder_formula(db_met)
    cls_h, cls_d = cofactor_class(host_met), cofactor_class(db_met)
    same_cofactor = bool(cls_h) and cls_h == cls_d

    if same_cofactor:
        # Redox partners share a name and a class; letting them merge would short the
        # couple and silently create a free electron source.
        state_h, state_d = _redox_state(host_met), _redox_state(db_met)
        if state_h and state_d and state_h != state_d:
            proposal.blocked = f"same cofactor but {state_h} vs {state_d}"
            add(Evidence("cofactor", f"redox states differ: {state_h} vs {state_d}",
                         _P_FORMULA_CONFLICT))
            proposal.score += _P_FORMULA_CONFLICT
        else:
            add(Evidence("cofactor",
                         f"both are {cls_h}, a protein-bound cofactor whose formula is a "
                         f"placeholder", _W_COFACTOR))
            proposal.score += _W_COFACTOR

    fh, fd = _formula_key(host_met), _formula_key(db_met)
    if fh and fd:
        if fh == fd:
            add(Evidence("formula",
                         f"same composition {host_met.formula}", _W_FORMULA))
            proposal.score += _W_FORMULA
        elif same_cofactor or placeholder:
            # Recorded, but neither rewarded nor penalised: the strings are placeholders
            # and the two databases were never going to agree on them.
            add(Evidence("formula",
                         f"formulas differ ({host_met.formula} vs {db_met.formula}) but "
                         f"are placeholders — no evidence either way", 0.0))
        else:
            add(Evidence("formula",
                         f"formulas differ: {host_met.formula} vs {db_met.formula}",
                         _P_FORMULA_CONFLICT))
            proposal.score += _P_FORMULA_CONFLICT
            if not proposal.blocked:
                proposal.blocked = "different molecular formulas"
    elif fh or fd:
        # The BiGG universal ships many metabolites without formulas. Absence is not a
        # mismatch; it is simply no evidence, and must not count against a pair.
        add(Evidence("formula", "one side has no formula — no evidence", 0.0))

    ch = getattr(host_met, "charge", None)
    cd = getattr(db_met, "charge", None)
    if ch is not None and cd is not None and ch != cd and not (same_cofactor or placeholder):
        add(Evidence("formula", f"charges differ: {ch} vs {cd}", _P_FORMULA_CONFLICT))
        proposal.score += _P_FORMULA_CONFLICT
        if not proposal.blocked:
            proposal.blocked = "different charges"

    # --- name -------------------------------------------------------------------------
    # Last resort, and equality only. Compared on both name and id so a database that
    # puts the human-readable string in the id (``beta_Carotene_c``) still matches.
    host_forms = {_plural_insensitive(normalise_name(x))
                  for x in (getattr(host_met, "name", ""), host_met.id) if x}
    db_forms = {_plural_insensitive(normalise_name(x))
                for x in (getattr(db_met, "name", ""), db_met.id) if x}
    host_forms.discard("")
    db_forms.discard("")
    if host_forms & db_forms:
        add(Evidence("name",
                     f"names agree once normalised: "
                     f"{sorted(host_forms & db_forms)[0]}", _W_NAME))
        proposal.score += _W_NAME

    return proposal


# --------------------------------------------------------------------------------------
# Generating candidates
# --------------------------------------------------------------------------------------
def _index_by(models: Iterable, key_fn) -> Dict[str, List]:
    index: Dict[str, List] = {}
    for met in models:
        key = key_fn(met)
        if key:
            index.setdefault(key, []).append(met)
    return index


def propose(host: cobra.Model, database: cobra.Model, *,
            only_unmatched: bool = True,
            metabolite_ids: Optional[Sequence[str]] = None,
            limit_per_host: int = 3,
            floor: float = REVIEW_FLOOR) -> List[Proposal]:
    """Suggest host↔database correspondences, best first.

    ``only_unmatched`` skips host metabolites the existing cross-reference machinery
    already resolves, which is the common case and not worth the user's attention. The
    point of this module is the residue those mechanisms cannot reach.

    ``metabolite_ids`` narrows the work to specific database metabolites — the natural
    call from a search result, where only the compounds a route actually uses matter.

    Neither model is modified.
    """
    db_mets = ([database.metabolites.get_by_id(m) for m in metabolite_ids
                if database.metabolites.has_id(m)]
               if metabolite_ids else list(database.metabolites))

    already: Set[str] = set()
    if only_unmatched:
        # Anything the token index already maps is out of scope here.
        try:
            already = set(namespace.match_database_to_host(host, database))
        except Exception:  # noqa: BLE001 — a failure here must not block reconciliation
            already = set()

    host_mets = list(host.metabolites)
    by_name = _index_by(host_mets, lambda m: _plural_insensitive(
        normalise_name(getattr(m, "name", "") or "")))
    by_id_name = _index_by(host_mets, lambda m: _plural_insensitive(
        normalise_name(m.id)))
    by_block = _index_by(host_mets, namespace._inchikey_block)
    by_key = _index_by(host_mets, _inchikey_full)
    by_token: Dict[str, List] = {}
    for met in host_mets:
        for tok in namespace.metabolite_tokens(met):
            if not tok.startswith("base:"):
                by_token.setdefault(tok, []).append(met)

    proposals: List[Proposal] = []
    for db_met in db_mets:
        if db_met.id in already:
            continue
        # Candidates from every channel, not only cross-references — the whole point,
        # since a disjoint annotation set yields no cross-reference candidate at all.
        candidates: Dict[str, object] = {}
        for source in (by_key.get(_inchikey_full(db_met)),
                       by_block.get(namespace._inchikey_block(db_met)),
                       by_name.get(_plural_insensitive(
                           normalise_name(getattr(db_met, "name", "") or ""))),
                       by_name.get(_plural_insensitive(normalise_name(db_met.id))),
                       by_id_name.get(_plural_insensitive(
                           normalise_name(getattr(db_met, "name", "") or ""))),
                       by_id_name.get(_plural_insensitive(normalise_name(db_met.id)))):
            for met in source or ():
                candidates[met.id] = met
        for tok in namespace.metabolite_tokens(db_met):
            if tok.startswith("base:"):
                continue
            for met in by_token.get(tok, ()):
                candidates[met.id] = met

        scored = [score_pair(met, db_met) for met in candidates.values()]
        scored = [p for p in scored if p.score >= floor or p.blocked]
        scored.sort(key=lambda p: (-p.score, p.host_id))
        proposals.extend(scored[:max(1, limit_per_host)])

    proposals.sort(key=lambda p: (-p.score, p.db_id, p.host_id))
    return proposals


def auto_accepted(proposals: Sequence[Proposal]) -> Dict[str, str]:
    """The subset safe to apply without review: ``{database_id: host_id}``."""
    mapping: Dict[str, str] = {}
    for p in proposals:
        if p.auto_acceptable and p.db_id not in mapping:
            mapping[p.db_id] = p.host_id
    return mapping


# --------------------------------------------------------------------------------------
# Persistence
# --------------------------------------------------------------------------------------
#: Bumped when the stored shape changes so an old project is recognised rather than
#: misread. A mapping that silently loses its evidence is worse than one that refuses to
#: load, because the numbers it produced would still look right.
MAPPING_SCHEMA = 1

#: Key under ``Project.datasets`` where the accepted mapping lives.
PROJECT_KEY = "id_reconciliation"


def to_record(accepted: Sequence[Proposal], *, database: str = "",
              overrides: Optional[Dict[str, str]] = None) -> dict:
    """Serialise accepted correspondences, keeping the reason each one was accepted.

    The evidence is stored, not just the pair. A mapping a reader cannot interrogate is
    a set of unexplained merges, and the whole point of the review step is that a human
    agreed to each one for a stated reason.
    """
    entries = []
    for p in accepted:
        entries.append({
            "database_id": p.db_id,
            "host_id": p.host_id,
            "score": round(float(p.score), 1),
            "confidence": p.confidence,
            "evidence": [{"kind": e.kind, "detail": e.detail, "weight": e.weight}
                         for e in p.evidence],
            "overridden": bool(p.blocked),
            "override_reason": (overrides or {}).get(p.db_id, ""),
            "blocked_reason": p.blocked,
        })
    entries.sort(key=lambda d: d["database_id"])
    return {"schema": MAPPING_SCHEMA, "database": database, "entries": entries}


def from_record(record: Optional[dict]) -> Dict[str, str]:
    """``{database_id: host_id}`` from a stored record, or empty if it is unusable."""
    if not isinstance(record, dict):
        return {}
    if int(record.get("schema", 0)) != MAPPING_SCHEMA:
        return {}
    mapping: Dict[str, str] = {}
    for entry in record.get("entries", []):
        db_id, host_id = entry.get("database_id"), entry.get("host_id")
        if db_id and host_id:
            mapping[str(db_id)] = str(host_id)
    return mapping


def describe_record(record: Optional[dict]) -> List[str]:
    """Human-readable lines for a report, one per accepted correspondence."""
    if not isinstance(record, dict) or not record.get("entries"):
        return []
    lines = []
    for entry in record["entries"]:
        why = "; ".join(e.get("detail", "") for e in entry.get("evidence", []))
        mark = " (user override)" if entry.get("overridden") else ""
        lines.append(f"{entry['database_id']} treated as {entry['host_id']}"
                     f"{mark}: {why}")
    return lines
