"""Read a metabolite's cross-references whatever spelling the database used.

Model files disagree about how to write an annotation key. BiGG-derived models use
``"KEGG Compound"``, ``"InChI Key"`` and ``"CHEBI"``, holding identifiers.org **URLs**;
the MIRIAM convention is ``kegg.compound``, ``inchikey``, ``chebi`` holding bare
accessions; others use ``kegg_id``, ``KEGG``, ``inchi_key``. Code that looks for one
spelling silently finds nothing in a model that uses another — and "nothing" is
indistinguishable from "this compound has no identifiers", so the failure is invisible.

That bug cost the thermodynamics analysis dearly: it resolved **zero** of the first 400
metabolites in the bundled database, so every compound fell through to a per-compound
online structure lookup. The analysis was slow because it was doing hundreds of web
requests, and scored almost nothing because most of those requests could not succeed.

One normalisation, used by everything that needs an identifier.
"""

from __future__ import annotations

import re
from typing import Dict, List

#: Canonical slug → the spellings that mean it, once punctuation and case are stripped.
_ALIASES = {
    "kegg": ("keggcompound", "kegg", "keggid", "keggligand", "keggc"),
    "chebi": ("chebi", "chebiid"),
    "inchikey": ("inchikey",),
    "inchi": ("inchi",),
    "smiles": ("smiles",),
    "metanetx": ("metanetxmnxchemical", "metanetxchemical", "metanetx", "mnx",
                 "mnxchemical"),
    "seed": ("seedcompound", "seed", "modelseed", "seedid"),
    "bigg": ("biggmetabolite", "bigg", "biggid"),
    "hmdb": ("humanmetabolomedatabase", "hmdb"),
    "metacyc": ("biocyc", "metacyc", "metacyccompound"),
    "lipidmaps": ("lipidmaps",),
    "reactome": ("reactome", "reactomecompound"),
}


def _slug(key: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(key).lower())


def _bare(value: str) -> str:
    """Strip an identifiers.org URL down to the accession it wraps."""
    text = str(value).strip()
    if "://" in text:
        text = text.rstrip("/").rsplit("/", 1)[-1]
    return text.strip()


def normalised(met) -> Dict[str, List[str]]:
    """``{canonical slug: [accession, …]}`` for a cobra metabolite.

    Every value is a bare accession — no URLs — and every alias spelling collapses onto
    the same slug, so callers can ask for ``"kegg"`` without knowing how the file wrote it.
    """
    ann = getattr(met, "annotation", None) or {}
    if not isinstance(ann, dict):
        return {}
    by_slug: Dict[str, List[str]] = {}
    for key, value in ann.items():
        slug = _slug(key)
        canonical = None
        for name, aliases in _ALIASES.items():
            if slug in aliases:
                canonical = name
                break
        if canonical is None:
            continue
        values = value if isinstance(value, (list, tuple)) else [value]
        for item in values:
            bare = _bare(item)
            if not bare:
                continue
            bucket = by_slug.setdefault(canonical, [])
            if bare not in bucket:
                bucket.append(bare)
    return by_slug


def first(met, slug: str) -> str:
    """The first accession of one kind, or ""."""
    found = normalised(met).get(slug) or []
    return found[0] if found else ""
