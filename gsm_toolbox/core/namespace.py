"""Cross-namespace metabolite matching.

Different databases identify the same compound with different ids: BiGG calls
glucose ``glc__D``, KEGG ``C00031``, MetaNetX ``MNXM41``. A heterologous-pathway
gap-fill can only connect a database reaction to the host network when they
share metabolite identifiers, so a KEGG/MetaNetX database fetched against a BiGG
host would otherwise find nothing.

This module matches metabolites *across* identifier systems by collecting every
cross-reference each metabolite carries (its own id, interpreted by namespace,
plus all the database links in its ``annotation``) and matching on any shared
reference, compartment-aware. It can then *translate* a database into the host's
namespace so the shared metabolites unify and gap-filling works.
"""

from __future__ import annotations

import re
from typing import Dict, Optional

import cobra

# annotation key (lower-cased) -> canonical namespace
_ANNOTATION_NS = {
    "kegg.compound": "kegg", "kegg": "kegg", "kegg.glycan": "kegg",
    "bigg.metabolite": "bigg", "bigg": "bigg",
    "metanetx.chemical": "metanetx", "metanetx": "metanetx", "mnx": "metanetx",
    "chebi": "chebi",
    "seed.compound": "seed", "seed": "seed",
    "inchi_key": "inchikey", "inchikey": "inchikey", "inchikey_first": "inchikey",
    "biocyc": "biocyc", "metacyc": "biocyc", "metacyc.compound": "biocyc",
    "hmdb": "hmdb",
    "sabiork": "sabiork", "sabiork.compound": "sabiork",
    "lipidmaps": "lipidmaps",
}


def _norm(ns: str, value: str) -> str:
    value = str(value).strip()
    if not value:
        return ""
    if ns == "chebi":
        v = value.upper().split("/")[-1]
        if not v.startswith("CHEBI:"):
            v = "CHEBI:" + v.split(":")[-1]
        return f"chebi:{v}"
    if ns == "inchikey":
        return f"inchikey:{value.upper()}"
    return f"{ns}:{value}"


def base_id(met: cobra.Metabolite) -> str:
    """Strip a metabolite's trailing compartment suffix (``glc__D_c`` -> ``glc__D``)."""
    mid = met.id
    comp = getattr(met, "compartment", "") or ""
    if comp and mid.endswith("_" + comp):
        return mid[: -(len(comp) + 1)]
    m = re.match(r"^(.*)_([a-z][a-z0-9]?)$", mid)
    return m.group(1) if m else mid


def metabolite_tokens(met: cobra.Metabolite) -> set:
    """All cross-reference tokens (``ns:value``) a metabolite carries."""
    tokens = set()
    ann = getattr(met, "annotation", None)
    if isinstance(ann, dict):
        for key, value in ann.items():
            ns = _ANNOTATION_NS.get(str(key).lower())
            if not ns:
                continue
            values = value if isinstance(value, (list, tuple, set)) else [value]
            for v in values:
                tok = _norm(ns, v)
                if tok:
                    tokens.add(tok)
    # Infer a token from the id itself, interpreted by its pattern.
    bid = base_id(met)
    if re.fullmatch(r"MNXM\d+", bid):
        tokens.add(_norm("metanetx", bid))
    elif re.fullmatch(r"[CG]\d{5}", bid):
        tokens.add(_norm("kegg", bid))
    elif re.fullmatch(r"cpd\d+", bid):
        tokens.add(_norm("seed", bid))
    elif bid:
        tokens.add(_norm("bigg", bid))
    return tokens


def _host_index(host: cobra.Model) -> Dict[str, Dict[str, str]]:
    """token -> {compartment: host_metabolite_id}."""
    index: Dict[str, Dict[str, str]] = {}
    for met in host.metabolites:
        comp = getattr(met, "compartment", "") or ""
        for tok in metabolite_tokens(met):
            index.setdefault(tok, {}).setdefault(comp, met.id)
    return index


def _formula_key(met) -> Optional[tuple]:
    """A comparable elemental composition, or None when the formula is missing."""
    try:
        el = getattr(met, "elements", None) or {}
        if not el:
            return None
        return tuple(sorted((str(k), int(v)) for k, v in el.items()))
    except Exception:  # noqa: BLE001 — an unparseable formula is simply 'unknown'
        return None


def _inchikey_block(met) -> Optional[str]:
    """The connectivity layer of a metabolite's InChIKey (skeleton identity)."""
    ann = getattr(met, "annotation", None) or {}
    for key in ("inchi_key", "inchikey", "InChI Key"):
        v = ann.get(key)
        if isinstance(v, (list, tuple)):
            v = v[0] if v else None
        if v:
            return str(v).strip().split("/")[-1].split("-", 1)[0].upper()
    return None


def _candidate_hosts(met, index: Dict[str, Dict[str, str]]) -> Dict[str, int]:
    """Host metabolites this database metabolite's cross-references point at, with the
    number of shared tokens supporting each."""
    comp = getattr(met, "compartment", "") or ""
    votes: Dict[str, int] = {}
    for tok in sorted(metabolite_tokens(met)):        # sorted ⇒ deterministic
        cands = index.get(tok)
        if not cands:
            continue
        hid = cands.get(comp) or cands.get("c") or sorted(cands.values())[0]
        if hid:
            votes[hid] = votes.get(hid, 0) + 1
    return votes


def _rank_candidates(met, votes: Dict[str, int], host_by_id) -> list:
    """Rank candidate host metabolites by CHEMICAL identity, then cross-reference support.

    Public databases contain contaminated entries: the merged universal's ``CO2``
    carries cobalt(2+) identifiers (KEGG C00175, ChEBI 48827/48828, HMDB00608 — BioCyc
    writes cobalt as ``CO+2``, which a naive merge fused with ``CO2``). Matching on the
    first cross-reference that happens to hit therefore mapped carbon dioxide onto
    cobalt, or raised a spurious ambiguity warning.

    Formula settles it: CO2 vs Co is decisive, and disagreement is strong evidence
    *against* a match even when identifiers agree. InChIKey skeleton is next, then the
    number of supporting cross-references. Returns ``[(score, host_id), …]`` best first.
    """
    fk, ik = _formula_key(met), _inchikey_block(met)
    scored = []
    for hid, n in votes.items():
        h = host_by_id.get(hid)
        score = float(n)                     # baseline: cross-reference support
        if h is not None:
            hfk = _formula_key(h)
            if fk and hfk:
                score += 100.0 if fk == hfk else -100.0
            hik = _inchikey_block(h)
            if ik and hik:
                score += 50.0 if ik == hik else -25.0
        scored.append((score, hid))
    scored.sort(key=lambda t: (-t[0], t[1]))          # deterministic tie-break by id
    return scored


def match_database_to_host(host: cobra.Model, database: cobra.Model) -> Dict[str, str]:
    """Map ``database`` metabolite ids -> host metabolite ids (compartment-aware).

    Only confident cross-reference matches are returned; unmatched database
    metabolites (e.g. a novel target and its intermediates) are omitted. When several
    host metabolites are implicated, the chemically consistent one wins (see
    :func:`_rank_candidates`) rather than whichever cross-reference was seen first.
    """
    index = _host_index(host)
    host_by_id = {m.id: m for m in host.metabolites}
    mapping: Dict[str, str] = {}
    for met in database.metabolites:
        votes = _candidate_hosts(met, index)
        if not votes:
            continue
        ranked = _rank_candidates(met, votes, host_by_id)
        if ranked:
            mapping[met.id] = ranked[0][1]
    return mapping


def count_shared(host: cobra.Model, database: cobra.Model) -> int:
    """How many database metabolites can be matched to the host."""
    return len(match_database_to_host(host, database))


def translation_report(host: cobra.Model, database: cobra.Model,
                       mapping: Optional[Dict[str, str]] = None,
                       *, relevant_ids: Optional[set] = None) -> Dict[str, object]:
    """Diagnostics for a database→host namespace match (Issue 14).

    Returns a dict with ``matched``/``unmatched`` counts and lists of
    ``ambiguous`` matches (a database metabolite whose cross-references point at
    more than one distinct host metabolite) and ``collisions`` (several database
    metabolites collapsing onto the same host metabolite). When ``relevant_ids``
    is given (e.g. the metabolites in a chosen route), the ambiguity/collision
    reports are restricted to those, so the user sees only what affects them.
    """
    if mapping is None:
        mapping = match_database_to_host(host, database)
    index = _host_index(host)

    host_by_id = {m.id: m for m in host.metabolites}
    ambiguous = []       # (db_id, [host_id, …])
    for met in database.metabolites:
        if relevant_ids is not None and met.id not in relevant_ids:
            continue
        votes = _candidate_hosts(met, index)
        if len(votes) < 2:
            continue
        # Several host metabolites are implicated — but that is only worth reporting if
        # chemistry cannot separate them. Contaminated public entries (CO2 carrying
        # cobalt identifiers) are resolved decisively by formula, and used to produce a
        # scary warning on almost every pathway. Warn only on a genuine tie.
        ranked = _rank_candidates(met, votes, host_by_id)
        if len(ranked) > 1 and ranked[0][0] - ranked[1][0] < 1e-9:
            ambiguous.append((met.id, sorted(h for _s, h in ranked if h)))

    # collisions: >1 database id mapped onto the same host id
    reverse: Dict[str, list] = {}
    for db_id, host_id in mapping.items():
        if relevant_ids is not None and db_id not in relevant_ids:
            continue
        reverse.setdefault(host_id, []).append(db_id)
    collisions = [(h, sorted(v)) for h, v in reverse.items() if len(v) > 1]

    total = (len(relevant_ids) if relevant_ids is not None else len(database.metabolites))
    matched = sum(1 for db_id in mapping
                  if relevant_ids is None or db_id in relevant_ids)
    return {
        "matched": matched,
        "unmatched": max(0, total - matched),
        "ambiguous": ambiguous,
        "collisions": collisions,
    }


def translate_database(host: cobra.Model, database: cobra.Model,
                       mapping: Optional[Dict[str, str]] = None) -> cobra.Model:
    """Return a copy of ``database`` with matched metabolites renamed to host ids.

    After translation, metabolites shared with the host have identical ids, so
    COBRApy gap-filling can connect the database's reactions to the host network.
    Unmatched metabolites keep their original id (they remain addable as novel
    species).
    """
    if mapping is None:
        mapping = match_database_to_host(host, database)

    new = cobra.Model(f"{database.id}_host_ns")
    mets: Dict[str, cobra.Metabolite] = {}

    def _target(met: cobra.Metabolite) -> cobra.Metabolite:
        tid = mapping.get(met.id, met.id)
        if tid not in mets:
            if host.metabolites.has_id(tid):
                hm = host.metabolites.get_by_id(tid)
                nm = cobra.Metabolite(tid, name=hm.name, formula=hm.formula,
                                      charge=hm.charge, compartment=hm.compartment)
            else:
                nm = cobra.Metabolite(tid, name=met.name, formula=met.formula,
                                      charge=met.charge, compartment=met.compartment)
            nm.annotation = dict(getattr(met, "annotation", {}) or {})
            mets[tid] = nm
        return mets[tid]

    reactions = []
    for rxn in database.reactions:
        nr = cobra.Reaction(rxn.id, name=rxn.name,
                            subsystem=getattr(rxn, "subsystem", ""))
        nr.bounds = rxn.bounds
        coeffs: Dict[cobra.Metabolite, float] = {}
        for met, c in rxn.metabolites.items():
            tgt = _target(met)
            coeffs[tgt] = coeffs.get(tgt, 0.0) + c
        coeffs = {m: c for m, c in coeffs.items() if c != 0.0}
        if not coeffs:
            continue
        nr.add_metabolites(coeffs)
        nr.annotation = dict(getattr(rxn, "annotation", {}) or {})
        reactions.append(nr)
    new.add_reactions(reactions)
    return new
