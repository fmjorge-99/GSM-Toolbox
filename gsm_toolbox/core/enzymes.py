"""EC-number suggestion and enzyme-candidate selection for a designed reaction.

Two related problems this module solves:

**1. Missing EC numbers.** Most reactions in the merged universal carry no EC annotation,
so the "EC" column is usually blank and the user has no starting point for finding an
enzyme. But the reaction almost always carries a *cross-reference* (KEGG, Rhea, MetaCyc,
SEED) from which an EC can be looked up, and when even that is absent, a reaction with the
same participants elsewhere in the loaded databases often does have one. `suggest_ec_numbers`
tries each of those in turn and says where every suggestion came from, so the user can
judge it rather than trust it blindly.

**2. "A reaction is needed" → "clone this gene" (Selenzyme-style selection).** The
reference tool for this is Selenzyme, which ranks candidate enzyme sequences by reaction
similarity. Selenzyme is an external web service and is not reachable from every network
(and is not bundled), so rather than depend on it, `enzyme_candidates` implements the
capability on open, reliable endpoints: EC → UniProtKB, ranked by review status
(SwissProt first) and by taxonomic proximity to the production host. That covers the
question a bench scientist actually asks — *which sequence do I clone?* — without a
single-point-of-failure dependency. See ``docs/enzyme_selection.md``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

_EC_RE = re.compile(r"\b\d+\.\d+\.\d+\.(?:\d+|-)\b")


@dataclass
class ECSuggestion:
    ec: str
    source: str          # where it came from, shown to the user
    confidence: str      # "annotated" | "cross-reference" | "inferred"
    detail: str = ""

    @property
    def url(self) -> str:
        """A canonical page for this EC number (ExPASy ENZYME)."""
        return f"https://enzyme.expasy.org/EC/{self.ec}"

    @property
    def brenda_url(self) -> str:
        return f"https://www.brenda-enzymes.org/enzyme.php?ecno={self.ec}"


@dataclass
class EnzymeCandidate:
    accession: str
    protein: str
    organism: str
    length: str = ""
    reviewed: str = ""
    score: float = 0.0
    why: str = ""

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.accession}/entry"


# --------------------------------------------------------------------------------------
# 1. EC suggestion
# --------------------------------------------------------------------------------------
def _annotated_ecs(rxn) -> List[str]:
    from .databases import reaction_ec_numbers
    try:
        return reaction_ec_numbers(rxn)
    except Exception:  # noqa: BLE001
        return []


def _xref(rxn, *keys) -> Optional[str]:
    ann = getattr(rxn, "annotation", None) or {}
    for k in keys:
        v = ann.get(k)
        if isinstance(v, (list, tuple)):
            v = v[0] if v else None
        if v:
            return str(v).strip().split("/")[-1]
    return None


def _ec_from_kegg(reaction_id: str, timeout: int = 20) -> List[str]:
    """KEGG's ENZYME line for a reaction id (``R05380`` → ``1.14.17.4``)."""
    from .databases import _http_get, KEGG_REST
    rid = reaction_id if reaction_id.startswith("R") else f"R{reaction_id}"
    try:
        text = _http_get(f"{KEGG_REST}/get/rn:{rid}", timeout=timeout)
    except Exception:  # noqa: BLE001
        return []
    out: List[str] = []
    for line in text.splitlines():
        if line.startswith("ENZYME"):
            out.extend(_EC_RE.findall(line))
    return out


def _ec_from_rhea(rhea_id: str, timeout: int = 25) -> List[str]:
    """Rhea's EC column for a RHEA id."""
    from .databases import _rhea_tsv
    q = rhea_id if rhea_id.upper().startswith("RHEA") else f"RHEA:{rhea_id}"
    try:
        rows = _rhea_tsv(q, columns="rhea-id,ec", limit=5, timeout=timeout)
    except Exception:  # noqa: BLE001
        return []
    out: List[str] = []
    for row in rows:
        if len(row) > 1:
            out.extend(_EC_RE.findall(row[1]))
    return out


def _reaction_signature(rxn) -> frozenset:
    """Participants of a reaction, compartment-agnostic, for similarity matching."""
    bases = set()
    for met in rxn.metabolites:
        mid = met.id
        bases.add(mid.rsplit("_", 1)[0].lower() if "_" in mid else mid.lower())
    return frozenset(bases)


def _ec_from_similar_reaction(rxn, db) -> List[tuple]:
    """EC numbers of database reactions with the SAME participant set.

    A reaction lacking an EC is often the same chemistry as a well-annotated one under a
    different id (the merged universal holds several namespaces side by side). Returns
    ``[(ec, source_reaction_id), …]``.
    """
    if db is None:
        return []
    sig = _reaction_signature(rxn)
    if len(sig) < 2:
        return []
    out: List[tuple] = []
    for other in db.reactions:
        if other.id == rxn.id:
            continue
        if _reaction_signature(other) != sig:
            continue
        for ec in _annotated_ecs(other):
            out.append((ec, other.id))
        if len(out) >= 8:
            break
    return out


def suggest_ec_numbers(rxn, db=None, *, online: bool = True) -> List[ECSuggestion]:
    """Every EC number that can be justified for ``rxn``, best-evidence first.

    Order of evidence: (1) an EC already annotated on the reaction; (2) an EC looked up
    from the reaction's own KEGG/Rhea cross-reference; (3) an EC borrowed from a database
    reaction with an identical participant set. Every suggestion records its source so
    the user can see *why* it was offered.
    """
    out: List[ECSuggestion] = []
    seen = set()

    def add(ec: str, source: str, confidence: str, detail: str = "") -> None:
        ec = ec.strip()
        if not ec or ec in seen or not _EC_RE.fullmatch(ec):
            return
        seen.add(ec)
        out.append(ECSuggestion(ec=ec, source=source, confidence=confidence,
                                detail=detail))

    for ec in _annotated_ecs(rxn):
        add(ec, "annotated on this reaction", "annotated")

    if online:
        kegg_id = _xref(rxn, "kegg.reaction", "kegg")
        if not kegg_id and re.fullmatch(r"R\d{5}", str(rxn.id)):
            kegg_id = rxn.id
        if kegg_id:
            for ec in _ec_from_kegg(kegg_id):
                add(ec, f"KEGG reaction {kegg_id}", "cross-reference",
                    "looked up from this reaction's KEGG cross-reference")
        rhea_id = _xref(rxn, "rhea")
        if rhea_id:
            for ec in _ec_from_rhea(rhea_id):
                add(ec, f"Rhea {rhea_id}", "cross-reference",
                    "looked up from this reaction's Rhea cross-reference")

    for ec, src in _ec_from_similar_reaction(rxn, db):
        add(ec, f"same chemistry as {src}", "inferred",
            "a database reaction with identical participants carries this EC — verify "
            "the transformation really is the same")
    return out


# --------------------------------------------------------------------------------------
# 2. Enzyme candidates (Selenzyme-style: reaction → sequences to clone)
# --------------------------------------------------------------------------------------
# Rough taxonomic affinity for common expression/production hosts. Used only to RANK
# candidates — never to exclude one — so an unusual but correct enzyme is still shown.
_HOST_HINTS: Dict[str, List[str]] = {
    "synechocystis": ["synechocystis", "synechococcus", "cyanobacteri", "anabaena",
                      "nostoc", "thermosynechococcus"],
    "escherichia": ["escherichia", "shigella", "salmonella", "enterobact", "klebsiella"],
    "saccharomyces": ["saccharomyces", "candida", "yeast", "kluyveromyces", "pichia"],
    "pseudomonas": ["pseudomonas", "azotobacter"],
    "bacillus": ["bacillus", "geobacillus", "paenibacillus"],
    "corynebacterium": ["corynebacterium", "mycobacterium", "rhodococcus"],
}


def _host_group(host_name: str) -> List[str]:
    low = (host_name or "").lower()
    for key, hints in _HOST_HINTS.items():
        if key in low or any(h in low for h in hints):
            return hints
    return []


def enzyme_candidates(ec: str, *, host_name: str = "", limit: int = 25
                      ) -> List[EnzymeCandidate]:
    """Candidate enzyme sequences for an EC number, ranked for use in ``host_name``.

    This is the practical substitute for a Selenzyme call (see the module docstring):
    UniProtKB is queried for the EC, then candidates are ranked by
    (a) **reviewed/SwissProt** status — a curated entry is far likelier to be the real
    activity — and (b) **taxonomic proximity to the production host**, since a
    cyanobacterial enzyme usually folds and functions better in a cyanobacterium.
    Ranking never removes a candidate.
    """
    from .databases import uniprot_enzymes_for_ec
    try:
        # Not reviewed-only: a curated hit is preferred by the ranking below, but for
        # unusual chemistry the only real enzyme may be an unreviewed entry, and hiding
        # it would defeat the purpose.
        df = uniprot_enzymes_for_ec(ec, reviewed_only=False)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"UniProt lookup failed for EC {ec}: {exc}") from exc
    hints = _host_group(host_name)
    out: List[EnzymeCandidate] = []
    for _i, row in df.iterrows():
        organism = str(row.get("Organism", ""))
        reviewed = str(row.get("Reviewed", ""))
        score, why = 0.0, []
        if reviewed.lower().startswith("review"):
            score += 10.0
            why.append("curated (SwissProt)")
        low = organism.lower()
        if hints and any(h in low for h in hints):
            score += 5.0
            why.append(f"related to the host ({host_name})")
        out.append(EnzymeCandidate(
            accession=str(row.get("Accession", "")), protein=str(row.get("Protein", "")),
            organism=organism, length=str(row.get("Length", "")), reviewed=reviewed,
            score=score, why=", ".join(why) or "matches the EC number"))
    out.sort(key=lambda c: (-c.score, c.accession))
    return out[:limit]
