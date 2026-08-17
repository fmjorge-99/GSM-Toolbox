"""Selenzyme-style enzyme selection: find enzymes by REACTION SIMILARITY.

The question this answers is the one the EC route cannot: *my designed reaction has no EC
number — which enzyme could still perform it?* That is the situation for every
RetroRules-generated step, and it is exactly what Selenzyme (Carbonell et al., SYNBIOCHEM)
was built for: fingerprint the query reaction, find the most similar characterised
reactions, and return the enzymes that catalyse them.

Why this is a native re-implementation rather than the upstream tool
-------------------------------------------------------------------
Upstream Selenzyme could not be installed into this application:

* it ships **only as a Docker service** (`synbiochem/selenzyme` → `start_server.sh`
  builds `continuumio/anaconda3:5.2.0` and runs a Flask app on port 5000);
* that image pins **2017-era dependencies** (``rdkit=2017.09.3.0``, Python 3.6,
  ``biopython=1.76``) which cannot coexist with this app's RDKit 2026 / Python 3.11;
* the algorithm repository it clones (``pablocarb/selenzy`` v1.0) carries **no licence**,
  so its code cannot be vendored here;
* the hosted service at ``selenzyme.synbiochem.co.uk`` refused connections when probed.

The *method*, however, is published and not encumbered, and every input it needs is open
data. This module implements it directly:

1. **Reference set** — Rhea reaction SMILES (``rhea-reaction-smiles.tsv``), the enzymes
   that catalyse them (``rhea2uniprot_sprot.tsv``, ~397k mappings) and their EC numbers
   (``rhea2ec.tsv``). Rhea is CC-BY 4.0 and downloaded on demand (a few MB).
2. **Fingerprint** — RDKit structural reaction fingerprints, the same family Selenzyme
   uses.
3. **Search** — Tanimoto similarity of the query reaction against the reference set.
4. **Ranking** — reaction similarity first, then curated-entry and host-taxonomy
   preference, mirroring Selenzyme's "reaction similarity × host compatibility" score.

Everything is local after the one-time download, so it works offline and cannot be broken
by a third-party service going down.
"""
from __future__ import annotations

import os
import pickle
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from . import cache

RHEA_TSV = "https://ftp.expasy.org/databases/rhea/tsv"
FILES = {
    "smiles": f"{RHEA_TSV}/rhea-reaction-smiles.tsv",
    "uniprot": f"{RHEA_TSV}/rhea2uniprot_sprot.tsv",
    "ec": f"{RHEA_TSV}/rhea2ec.tsv",
}
INDEX_VERSION = 2          # bump when the stored index layout changes
# v2: difference fingerprints on cofactor-stripped reactions (v1 could not discriminate)

# Ubiquitous cofactors and small species. They are stripped from BOTH the reference
# reactions and the query before fingerprinting, for two reasons:
#   * a bulky cofactor (NAD, CoA, SAM) dominates a whole-molecule fingerprint, so two
#     completely different transformations that share NAD look similar;
#   * rule-generated (RetroRules) steps usually omit cofactors entirely, so a query would
#     otherwise never match the fully-written database reaction. Measured: ethanol →
#     acetaldehyde scored 0.089 against the real alcohol-dehydrogenase reaction with
#     cofactors left in — far below unrelated noise.
# What remains is the transformation the enzyme actually performs.
_COFACTOR_SMILES = [
    "O", "[H+]", "O=O", "O=C=O", "[NH4+]", "OO", "[OH-]",
    "OP(=O)(O)O", "OP(=O)(O)OP(=O)(O)O",            # phosphate, diphosphate
    "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",  # NAD+
    "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc45)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",   # NADH
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",  # ATP
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",           # ADP
    "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O",                    # AMP
    "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](OP(=O)(O)O)[C@@H]1O)[C@@H](O)C(=O)NCCC(=O)NCCS",  # CoA
    "C[S+](CC[C@H](N)C(=O)[O-])C[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",  # SAM
    "OC(=O)[C@@H](N)CCSC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",          # SAH
]
_COFACTOR_BLOCKS: Optional[set] = None


@dataclass
class SelenzymeHit:
    """One candidate enzyme, with the evidence that produced it."""

    accession: str
    similarity: float                 # Tanimoto vs the query reaction (0-1)
    rhea_id: str = ""
    ec: str = ""
    organism: str = ""
    protein: str = ""
    score: float = 0.0                # similarity + host/curation bonuses
    why: str = ""

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.accession}/entry"

    @property
    def rhea_url(self) -> str:
        rid = self.rhea_id.replace("RHEA:", "")
        return f"https://www.rhea-db.org/rhea/{rid}"


def _cofactor_blocks() -> set:
    """InChIKey connectivity blocks of the cofactors to strip (computed once)."""
    global _COFACTOR_BLOCKS
    if _COFACTOR_BLOCKS is None:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        blocks = set()
        for smi in _COFACTOR_SMILES:
            try:
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    blocks.add(Chem.MolToInchiKey(m).split("-", 1)[0])
            except Exception:  # noqa: BLE001
                pass
        _COFACTOR_BLOCKS = blocks
    return _COFACTOR_BLOCKS


def strip_cofactors(reaction_smiles: str) -> str:
    """Remove ubiquitous cofactors from both sides, leaving the core transformation.

    Returns "" when a side becomes empty (a reaction that is *only* cofactor chemistry
    carries no transformation worth matching).
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    if ">>" not in reaction_smiles:
        return ""
    drop = _cofactor_blocks()

    def keep(side: str) -> List[str]:
        out = []
        for part in side.split("."):
            part = part.strip()
            if not part:
                continue
            try:
                m = Chem.MolFromSmiles(part)
                if m is None:
                    continue
                if m.GetNumHeavyAtoms() <= 1:      # lone ions/atoms carry no chemistry
                    continue
                if Chem.MolToInchiKey(m).split("-", 1)[0] in drop:
                    continue
                out.append(Chem.MolToSmiles(m))
            except Exception:  # noqa: BLE001
                continue
        return out

    left_s, right_s = reaction_smiles.split(">>", 1)
    left, right = keep(left_s), keep(right_s)
    if not left or not right:
        return ""
    return ".".join(left) + ">>" + ".".join(right)


def _reaction_fp(reaction_smiles: str):
    """Difference fingerprint of the cofactor-stripped reaction, or None.

    The *difference* fingerprint encodes what the reaction CHANGES rather than what it
    contains, which is what makes two alcohol dehydrogenases look alike and an unrelated
    reaction look different (measured: a decoy scored 0.172 with a structural
    fingerprint but 0.010 with this one).
    """
    from rdkit.Chem import rdChemReactions
    core = strip_cofactors(reaction_smiles)
    if not core:
        return None
    try:
        rxn = rdChemReactions.ReactionFromSmarts(core, useSmiles=True)
        if rxn is None:
            return None
        return rdChemReactions.CreateDifferenceFingerprintForReaction(rxn)
    except Exception:  # noqa: BLE001
        return None


def data_dir() -> str:
    d = os.path.join(cache.base_dir(), "selenzyme")
    os.makedirs(d, exist_ok=True)
    return d


def _index_path() -> str:
    return os.path.join(data_dir(), f"rhea_rxnfp_v{INDEX_VERSION}.pkl")


def is_installed() -> bool:
    """True once the reference data has been downloaded and indexed."""
    p = _index_path()
    return os.path.exists(p) and os.path.getsize(p) > 10_000


def rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except Exception:  # noqa: BLE001
        return False


# --------------------------------------------------------------------------------------
# Install: download the open Rhea tables and build the fingerprint index
# --------------------------------------------------------------------------------------
def _download(url: str, dest: str, progress: Optional[Callable] = None) -> str:
    from urllib.request import urlopen, Request
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return dest
    if progress:
        progress(f"Downloading {os.path.basename(dest)}…")
    req = Request(url, headers={"User-Agent": "gsm-toolbox"})
    with urlopen(req, timeout=120) as fh, open(dest, "wb") as out:
        out.write(fh.read())
    return dest


def _parse_map(path: str, key_col: int, val_col: int, *, header: bool) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if header and i == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(key_col, val_col):
                continue
            out.setdefault(parts[key_col].strip(), []).append(parts[val_col].strip())
    return out


def _ec_for(rid: str, rhea2ec: Dict[str, List[str]]) -> str:
    """EC numbers for a Rhea reaction id.

    Rhea allocates ids in blocks of four — master N, then N+1 (L→R), N+2 (R→L),
    N+3 (bidirectional). ``rhea-reaction-smiles.tsv`` uses the *directional* ids while
    ``rhea2ec.tsv`` is keyed on the *master*, so a direct join silently finds nothing.
    Fall back through the block to the master.
    """
    for cand in (rid, *(str(int(rid) - k) for k in (1, 2, 3) if rid.isdigit())):
        vals = rhea2ec.get(cand)
        if vals:
            return ",".join(sorted(set(vals)))
    return ""


def install(progress: Optional[Callable] = None, *, limit: Optional[int] = None) -> dict:
    """Download the Rhea reference data and build the reaction-fingerprint index.

    One-time, a few MB, and everything afterwards runs locally. ``progress`` is an
    optional ``callable(message)``. Returns a summary dict.
    """
    if not rdkit_available():
        raise RuntimeError("RDKit is required for reaction-similarity search.")
    from rdkit import RDLogger
    from rdkit.Chem import rdChemReactions
    RDLogger.DisableLog("rdApp.*")          # Rhea SMILES include exotic species

    d = data_dir()
    paths = {k: _download(v, os.path.join(d, os.path.basename(v)), progress)
             for k, v in FILES.items()}

    if progress:
        progress("Reading enzyme and EC mappings…")
    # rhea2uniprot_sprot / rhea2ec: RHEA_ID, DIRECTION, MASTER_ID, ID  (with header)
    rhea2up = _parse_map(paths["uniprot"], 0, 3, header=True)
    rhea2ec = _parse_map(paths["ec"], 0, 3, header=True)

    if progress:
        progress("Building reaction fingerprints (one-off, a few minutes)…")
    fps, meta = [], []
    n_read = n_kept = 0
    with open(paths["smiles"], "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rid, smi = parts[0].strip(), parts[1].strip()
            n_read += 1
            # Only reactions with a known enzyme are useful as an answer.
            accs = rhea2up.get(rid)
            if not accs:
                continue
            fp = _reaction_fp(smi)
            if fp is None:      # unparseable, or pure-cofactor chemistry
                continue
            fps.append(fp)
            meta.append((rid, tuple(sorted(set(accs))[:40]), _ec_for(rid, rhea2ec)))
            n_kept += 1
            if limit and n_kept >= limit:
                break
            if progress and n_kept % 2000 == 0:
                progress(f"Fingerprinted {n_kept} reactions…")

    if not fps:
        raise RuntimeError("No Rhea reactions could be fingerprinted.")
    if progress:
        progress("Saving index…")
    with open(_index_path(), "wb") as out:
        pickle.dump({"version": INDEX_VERSION, "fps": fps, "meta": meta}, out,
                    protocol=pickle.HIGHEST_PROTOCOL)
    return {"reactions_read": n_read, "reactions_indexed": n_kept,
            "enzyme_mappings": sum(len(v) for v in rhea2up.values()),
            "index_path": _index_path(),
            "index_mb": round(os.path.getsize(_index_path()) / 1e6, 1)}


_INDEX = None


def _load_index():
    global _INDEX
    if _INDEX is None:
        if not is_installed():
            raise RuntimeError(
                "The reaction-similarity reference data is not installed. Enable it in "
                "Settings ▸ Preferences ▸ Enable Selenzyme (reaction-similarity enzyme "
                "search).")
        with open(_index_path(), "rb") as fh:
            _INDEX = pickle.load(fh)
    return _INDEX


# --------------------------------------------------------------------------------------
# Query
# --------------------------------------------------------------------------------------
def reaction_smiles_for(rxn) -> str:
    """Build ``reactants>>products`` SMILES for a cobra reaction, or "" if structures
    are unavailable (most genome-scale metabolites carry only ids, not structures)."""
    from .retrorules import metabolite_smiles
    left, right = [], []
    for met, coeff in rxn.metabolites.items():
        smi = ""
        try:
            smi = metabolite_smiles(met) or ""
        except Exception:  # noqa: BLE001
            smi = ""
        if not smi:
            return ""                     # a partial reaction would be misleading
        n = max(1, int(abs(coeff)))
        (left if coeff < 0 else right).extend([smi] * n)
    if not left or not right:
        return ""
    return ".".join(left) + ">>" + ".".join(right)


def find_enzymes(reaction_smiles: str, *, host_name: str = "", limit: int = 25,
                 min_similarity: float = 0.1) -> List[SelenzymeHit]:
    """Enzymes for a reaction, ranked by reaction similarity then host compatibility.

    ``reaction_smiles`` is ``reactants>>products``. Returns at most ``limit`` hits with
    similarity ≥ ``min_similarity``; an empty list means nothing in Rhea resembles this
    chemistry closely enough to suggest an enzyme honestly.
    """
    if not rdkit_available():
        raise RuntimeError("RDKit is required for reaction-similarity search.")
    from rdkit import DataStructs, RDLogger
    from rdkit.Chem import rdChemReactions
    RDLogger.DisableLog("rdApp.*")

    idx = _load_index()
    qfp = _reaction_fp(reaction_smiles)
    if qfp is None:
        raise RuntimeError(
            "This reaction could not be fingerprinted — after removing cofactors there "
            "is no transformation left to compare, or its structures are unparseable.")

    # Difference fingerprints are sparse int vectors; Tanimoto is applied pairwise.
    ref = idx["fps"]
    sims = [DataStructs.TanimotoSimilarity(qfp, f) for f in ref]
    order = sorted(range(len(sims)), key=lambda i: -sims[i])

    from .enzymes import _host_group
    hints = _host_group(host_name)
    hits: List[SelenzymeHit] = []
    seen_acc = set()
    for i in order:
        s = float(sims[i])
        if s < min_similarity:
            break
        rid, accs, ec = idx["meta"][i]
        for acc in accs:
            if acc in seen_acc:
                continue
            seen_acc.add(acc)
            why = [f"reaction similarity {s:.2f}"]
            hits.append(SelenzymeHit(
                accession=acc, similarity=round(s, 4), rhea_id=f"RHEA:{rid}", ec=ec,
                score=s, why="; ".join(why)))
            if len(hits) >= limit * 3:      # gather, then annotate/rank a shortlist
                break
        if len(hits) >= limit * 3:
            break

    # Annotate the shortlist with organism/protein from UniProt and apply host preference.
    _annotate(hits[:limit * 2], hints, host_name)
    hits.sort(key=lambda h: (-h.score, -h.similarity, h.accession))
    return hits[:limit]


def _annotate(hits: List[SelenzymeHit], hints: List[str], host_name: str) -> None:
    """Fill in organism/protein for the shortlist and add the host-compatibility bonus.

    One batched UniProt query keeps this to a single round-trip.
    """
    if not hits:
        return
    accs = [h.accession for h in hits]
    info: Dict[str, Tuple[str, str, str]] = {}
    try:
        import urllib.parse
        from .databases import _http_get, UNIPROT_REST
        query = " OR ".join(f"accession:{a}" for a in accs[:100])
        params = urllib.parse.urlencode({
            "query": query, "format": "tsv", "size": min(100, len(accs)),
            "fields": "accession,protein_name,organism_name,reviewed"})
        text = _http_get(f"{UNIPROT_REST}?{params}", timeout=60)
        for ln in text.splitlines()[1:]:
            p = ln.split("\t")
            if len(p) >= 4:
                info[p[0].strip()] = (p[1], p[2], p[3])
    except Exception:  # noqa: BLE001 — annotation is a nicety; similarity already stands
        pass
    for h in hits:
        prot, org, rev = info.get(h.accession, ("", "", ""))
        h.protein, h.organism = prot, org
        extra = []
        if rev.lower().startswith("review"):
            h.score += 0.15
            extra.append("curated (SwissProt)")
        if hints and org and any(x in org.lower() for x in hints):
            h.score += 0.25
            extra.append(f"related to the host ({host_name})")
        if extra:
            h.why += "; " + "; ".join(extra)
