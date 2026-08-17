"""Access to large external reaction/enzyme databases for Pathway Design.

Three online resources are wrapped here, all free for academic use:

* **MetaNetX / MNXref** — a reconciled universal reaction database (KEGG +
  MetaCyc + BiGG + SEED + Rhea). We download the compact ``reac_prop.tsv``
  (~10 MB; the compound files are 600-800 MB and are intentionally avoided) and
  build a COBRA "universal" model. EC numbers ship in the file. MetaCyc content
  reaches the user through this reconciliation (MetaCyc itself has no free bulk
  API).
* **KEGG** — queried live through its REST API for *targeted* discovery: resolve
  a product name (e.g. "isoprene") to a compound, pull the reactions that make
  or consume it, and build a focused universal of a few hundred reactions with
  names and EC numbers. Good when you know the target by name.
* **UniProt** — given an EC number (from a predicted reaction) list candidate
  enzymes, their source organisms and sequence links, i.e. genes you could
  express to realise that step.

Everything here is pure Python (no Qt) and network access is wrapped so the GUI
can surface friendly errors and run downloads on a worker thread.
"""

from __future__ import annotations

import os
import re
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import cobra
import pandas as pd

KEGG_REST = "https://rest.kegg.jp"
RHEA_REST = "https://www.rhea-db.org/rhea"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/search"
METANETX_REAC_URL = "https://www.metanetx.org/ftp/latest/reac_prop.tsv"
METANETX_CHEM_XREF_URL = "https://www.metanetx.org/ftp/latest/chem_xref.tsv"
METANETX_REAC_XREF_URL = "https://www.metanetx.org/ftp/latest/reac_xref.tsv"

_MODELSEED_BASE = ("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/"
                   "master/Biochemistry")
MODELSEED_REACTIONS_URL = f"{_MODELSEED_BASE}/reactions.tsv"
MODELSEED_COMPOUNDS_URL = f"{_MODELSEED_BASE}/compounds.tsv"

# reac_xref source-prefix -> annotation key; the short reaction id is preferred
# from these sources in order (BiGG ids like "PAL" make the best model ids).
_REAC_XREF_SOURCE_NS = {
    "bigg.reaction": "bigg.reaction", "biggr": "bigg.reaction",
    "kegg.reaction": "kegg.reaction", "keggr": "kegg.reaction",
    "seed.reaction": "seed.reaction", "seedr": "seed.reaction",
    "metacyc.reaction": "metacyc.reaction", "metacycr": "metacyc.reaction",
    "rhea": "rhea", "rh": "rhea",
    "sabiork.reaction": "sabiork.reaction", "sabiorkr": "sabiork.reaction",
}
_REAC_SHORT_PRIORITY = ("bigg.reaction", "kegg.reaction", "seed.reaction", "metacyc.reaction")
_REAC_NAME_PRIORITY = ("kegg.reaction", "bigg.reaction", "metacyc.reaction")

# MetaNetX chem_xref source-prefix -> annotation key understood by the matcher
# (gsm_toolbox.core.namespace). Both the canonical and the "*M" model variants map
# to the same namespace.
_XREF_SOURCE_NS = {
    "chebi": "chebi",
    "kegg.compound": "kegg.compound", "keggc": "kegg.compound",
    "kegg.drug": "kegg.compound", "keggd": "kegg.compound",
    "bigg.metabolite": "bigg.metabolite", "biggm": "bigg.metabolite",
    "seed.compound": "seed.compound", "seedm": "seed.compound",
    "metacyc.compound": "metacyc.compound", "metacycm": "metacyc.compound",
    "hmdb": "hmdb",
    "sabiork.compound": "sabiork.compound", "sabiorkm": "sabiork.compound",
    "lipidmaps": "lipidmaps", "lipidmapsm": "lipidmaps",
}

# Annotation keys under which models commonly record EC numbers.
_EC_KEYS = ("ec-code", "ec_code", "ecnumbers", "ec number", "ec", "enzyme",
            "ec-number", "ec_number")
_EC_RE = re.compile(r"\b\d+\.\d+\.\d+\.\d+\b")


class DatabaseError(Exception):
    """Raised when an external database cannot be queried or parsed."""


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #
def _cache_dir() -> str:
    base = os.path.join(os.path.expanduser("~"), ".gsm_toolbox", "databases")
    os.makedirs(base, exist_ok=True)
    return base


def _cached_model_path(name: str) -> str:
    """Path under the persistent database cache for a built universal model."""
    return os.path.join(_cache_dir(), name)


def _save_built_model(model: cobra.Model, name: str) -> None:
    """Persist a built universal model as JSON for instant reuse next time."""
    try:
        from cobra.io import save_json_model
        save_json_model(model, _cached_model_path(name))
    except Exception:  # noqa: BLE001 - caching is best-effort, never fatal
        pass


def _load_built_model(name: str) -> Optional[cobra.Model]:
    path = _cached_model_path(name)
    if not os.path.exists(path) or os.path.getsize(path) < 100:
        return None
    try:
        from .io_models import _load_json_tolerant
        return _load_json_tolerant(path)
    except Exception:  # noqa: BLE001 - fall back to rebuilding
        return None


def _http_get(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "GSM-ToolBox"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", "replace")
    except Exception as exc:  # noqa: BLE001
        raise DatabaseError(f"Request failed:\n{url}\n\n{exc}") from exc


def _download(url: str, dest: str, timeout: int = 300) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "GSM-ToolBox"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as fh:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    except Exception as exc:  # noqa: BLE001
        raise DatabaseError(f"Download failed:\n{url}\n\n{exc}") from exc
    return dest


def reaction_ec_numbers(rxn: cobra.Reaction) -> List[str]:
    """Best-effort extraction of EC numbers from a reaction's annotation."""
    found: List[str] = []
    ann = getattr(rxn, "annotation", None)
    if isinstance(ann, dict):
        for key, value in ann.items():
            if key.lower() not in _EC_KEYS:
                continue
            values = value if isinstance(value, (list, tuple, set)) else [value]
            for v in values:
                found.extend(_EC_RE.findall(str(v)))
    # de-duplicate, keep order
    seen = set()
    out = []
    for ec in found:
        if ec not in seen:
            seen.add(ec)
            out.append(ec)
    return out


# --------------------------------------------------------------------------- #
# KEGG (live, targeted)
# --------------------------------------------------------------------------- #
def _kegg(path: str) -> str:
    return _http_get(f"{KEGG_REST}/{path}")


def kegg_find_compound(name: str) -> List[Tuple[str, str]]:
    """Resolve a compound *name* to KEGG compound ids: ``[(C-number, names), …]``."""
    name = name.strip()
    if not name:
        return []
    text = _kegg(f"find/compound/{urllib.parse.quote(name)}")
    out: List[Tuple[str, str]] = []
    for line in text.splitlines():
        if "\t" not in line:
            continue
        cid, desc = line.split("\t", 1)
        out.append((cid.replace("cpd:", ""), desc.strip()))
    return out


def _kegg_list_names(ids: List[str], prefix: str = "cpd",
                     deadline: Optional[float] = None) -> Dict[str, str]:
    """Map KEGG ids -> primary name using batched ``/list`` calls (10 per call)."""
    names: Dict[str, str] = {}
    ids = [i for i in ids if i]
    for i in range(0, len(ids), 10):
        if deadline is not None and time.monotonic() > deadline:
            break
        batch = "+".join(ids[i:i + 10])
        try:
            text = _kegg(f"list/{batch}")
        except DatabaseError:
            continue
        for line in text.splitlines():
            if "\t" not in line:
                continue
            key, desc = line.split("\t", 1)
            key = key.split(":")[-1]
            names[key] = desc.split(";")[0].strip()
    return names


def kegg_reactions_for_compound(cid: str) -> List[str]:
    """KEGG reaction ids that involve a compound."""
    cid = cid.replace("cpd:", "").strip()
    text = _kegg(f"link/reaction/cpd:{cid}")
    rids = []
    for line in text.splitlines():
        if "\t" in line:
            rids.append(line.split("\t", 1)[1].replace("rn:", "").strip())
    return rids


@dataclass
class _KeggReaction:
    rid: str
    name: str
    equation: str          # raw KEGG equation (C-ids)
    ec: List[str] = field(default_factory=list)
    reversible: bool = True
    substrates: Dict[str, float] = field(default_factory=dict)
    products: Dict[str, float] = field(default_factory=dict)


def _parse_kegg_side(side: str) -> Dict[str, float]:
    coeffs: Dict[str, float] = {}
    for term in side.split(" + "):
        term = term.strip()
        if not term:
            continue
        parts = term.split()
        coef = 1.0
        cid = parts[-1]
        if len(parts) >= 2:
            try:
                coef = float(parts[0])
            except ValueError:
                coef = 1.0  # polymer counts like "(n)" -> treat as 1
        coeffs[cid] = coeffs.get(cid, 0.0) + coef
    return coeffs


def _parse_kegg_entry(text: str) -> Optional[_KeggReaction]:
    rid = name = equation = ""
    ec: List[str] = []
    for raw in text.splitlines():
        if raw.startswith("ENTRY"):
            rid = raw.split()[1]
        elif raw.startswith("NAME"):
            name = raw[len("NAME"):].strip()
        elif raw.startswith("EQUATION"):
            equation = raw[len("EQUATION"):].strip()
        elif raw.startswith("ENZYME"):
            ec = _EC_RE.findall(raw)
    if not rid or not equation:
        return None
    reversible = "<=>" in equation
    sep = "<=>" if reversible else ("=>" if "=>" in equation else "<=")
    lhs, _, rhs = equation.partition(sep)
    subs = _parse_kegg_side(lhs)
    prods = _parse_kegg_side(rhs)
    return _KeggReaction(rid=rid, name=name, equation=equation, ec=ec,
                         reversible=reversible, substrates=subs, products=prods)


def kegg_reactions_for_ec(ec: str) -> List[str]:
    """KEGG reaction ids catalysed by an EC number (accepts ``4.1.1.102`` or ``EC:4.1.1.102``)."""
    ec = ec.strip()
    for pre in ("EC:", "ec:", "EC ", "ec "):
        if ec.startswith(pre):
            ec = ec[len(pre):].strip()
    if not ec:
        return []
    text = _kegg(f"link/reaction/ec:{ec}")
    rids = []
    for line in text.splitlines():
        if "\t" in line:
            rids.append(line.split("\t", 1)[1].replace("rn:", "").strip())
    return rids


def kegg_compound_names(cids) -> Dict[str, str]:
    """Map KEGG compound ids -> primary name (batched)."""
    return _kegg_list_names([c.replace("cpd:", "") for c in cids])


def kegg_get_reactions(rids: List[str], deadline: Optional[float] = None) -> List[_KeggReaction]:
    """Fetch and parse KEGG reaction entries (batched 10 per request)."""
    out: List[_KeggReaction] = []
    rids = list(dict.fromkeys(r for r in rids if r))
    for i in range(0, len(rids), 10):
        if deadline is not None and time.monotonic() > deadline:
            break
        batch = "+".join(f"rn:{r}" for r in rids[i:i + 10])
        try:
            text = _kegg(f"get/{batch}")
        except DatabaseError:
            continue
        for entry in text.split("\n///"):
            entry = entry.strip()
            if not entry:
                continue
            parsed = _parse_kegg_entry(entry)
            if parsed is not None:
                out.append(parsed)
    return out


def build_kegg_pathway_db(
    target: str,
    *,
    expand_steps: int = 1,
    max_reactions: int = 300,
    time_budget: float = 45.0,
    max_frontier: int = 80,
) -> Tuple[cobra.Model, str, str]:
    """Build a focused universal model around a target product, from KEGG.

    ``target`` is a compound name (e.g. "isoprene") or a KEGG compound id
    (``C16521``). Reactions touching the target are collected; ``expand_steps``
    extra rounds pull in reactions of the neighbouring metabolites. To stay
    responsive on highly-connected targets (e.g. astaxanthin), expansion is
    bounded by ``max_reactions``, a per-round frontier cap (``max_frontier``) and
    a wall-clock ``time_budget`` — whatever has been gathered when a limit is hit
    is returned rather than hanging.
    """
    target = target.strip()
    if re.fullmatch(r"C\d{5}", target):
        seed_cid, label = target, target
    else:
        hits = kegg_find_compound(target)
        if not hits:
            raise DatabaseError(f"KEGG has no compound matching '{target}'.")
        # Prefer a hit whose synonyms contain an exact (case-insensitive) match
        # to the query, so "isoprene" picks C16521 over "cis-1,4-Polyisoprene".
        want = target.lower()
        best = None
        for cid, desc in hits:
            names = [n.strip().lower() for n in desc.split(";")]
            if want in names:
                best = (cid, desc)
                break
        cid, desc = best or hits[0]
        seed_cid = cid
        label = f"{desc.split(';')[0].strip()} ({seed_cid})"

    # Reuse a previously built KEGG database for this seed/neighbourhood.
    cache_name = f"kegg_{seed_cid}_s{expand_steps}_m{max_reactions}.json"
    cache_path = _cached_model_path(cache_name)
    cached = _load_built_model(cache_name)
    if cached is not None:
        return cached, f"KEGG: {label}", cache_path

    start = time.monotonic()
    deadline = start + time_budget

    def _over_budget() -> bool:
        return time.monotonic() > deadline

    collected: Dict[str, _KeggReaction] = {}
    frontier = {seed_cid}
    seen_compounds: set = set()
    for _ in range(max(1, expand_steps + 1)):
        rids: List[str] = []
        # Cap how many compounds we expand per round so a hub metabolite can't
        # trigger thousands of KEGG requests.
        for cid in sorted(frontier)[:max_frontier]:
            if cid in seen_compounds:
                continue
            seen_compounds.add(cid)
            if _over_budget():
                break
            try:
                rids.extend(kegg_reactions_for_compound(cid))
            except DatabaseError:
                continue
        rids = [r for r in dict.fromkeys(rids) if r not in collected][:max_reactions]
        if not rids:
            break
        for rxn in kegg_get_reactions(rids, deadline=deadline):
            if len(collected) >= max_reactions:
                break
            collected[rxn.rid] = rxn
        if len(collected) >= max_reactions or _over_budget():
            break
        frontier = {c for rxn in collected.values()
                    for c in (list(rxn.substrates) + list(rxn.products))}

    if not collected:
        raise DatabaseError(
            f"No KEGG reactions found around '{target}'. It may be too highly connected "
            "to explore, or not linked to reactions in KEGG — try a specific KEGG id or a "
            "different database.")

    # Resolve compound names for nicer metabolite labels.
    all_cids = {c for rxn in collected.values()
                for c in (list(rxn.substrates) + list(rxn.products))}
    # Name the SEED compound first, then the rest. Naming is batched and budgeted, and
    # when the budget runs out the remaining compounds keep their bare KEGG id as their
    # name. Alphabetical order therefore decided who got a name: fetching "erythromycin"
    # named C00001/C00007 but ran out before C01912 — so the database came back without
    # the very compound the user asked for being findable by name. The seed must never
    # lose that race, and the budget now scales with the amount to look up (this is a
    # one-off fetch that already takes a minute; 15s of naming was a false economy).
    ordered = [seed_cid] + sorted(all_cids - {seed_cid})
    budget = max(30.0, 0.2 * len(ordered))
    names = _kegg_list_names(ordered, deadline=time.monotonic() + budget)

    model = cobra.Model("kegg_pathway_db")
    mets: Dict[str, cobra.Metabolite] = {}

    def _met(cid: str) -> cobra.Metabolite:
        if cid not in mets:
            m = cobra.Metabolite(f"{cid}_c", name=names.get(cid, cid), compartment="c")
            # Keep the KEGG id as a cross-reference so the metabolite can still be
            # matched to a host model (and re-named to a readable id) later.
            m.annotation["kegg.compound"] = cid
            mets[cid] = m
        return mets[cid]

    reactions = []
    for rxn in collected.values():
        r = cobra.Reaction(rxn.rid, name=rxn.name or rxn.rid)
        r.lower_bound = -1000.0 if rxn.reversible else 0.0
        r.upper_bound = 1000.0
        r.annotation["kegg.reaction"] = rxn.rid
        coeffs = {}
        for cid, c in rxn.substrates.items():
            coeffs[_met(cid)] = -abs(c)
        for cid, c in rxn.products.items():
            coeffs[_met(cid)] = coeffs.get(_met(cid), 0.0) + abs(c)
        if not coeffs:
            continue
        r.add_metabolites(coeffs)
        if rxn.ec:
            r.annotation["ec-code"] = rxn.ec
        reactions.append(r)
    model.add_reactions(reactions)
    _save_built_model(model, cache_name)
    return model, f"KEGG: {label}", cache_path


# --------------------------------------------------------------------------- #
# Rhea (open-licence expert-curated reactions, ChEBI-based)
# --------------------------------------------------------------------------- #
# Currency/cofactor names skipped when expanding the compound frontier, so a hub
# metabolite (water, NAD…) can't trigger a combinatorial fan-out of queries.
_RHEA_CURRENCY = {
    "h2o", "water", "h(+)", "h+", "proton", "co2", "carbon dioxide", "o2", "dioxygen",
    "nad(+)", "nadh", "nadp(+)", "nadph", "atp", "adp", "amp", "diphosphate", "phosphate",
    "coa", "coenzyme a", "nh4(+)", "hydrogen peroxide", "fad", "fadh2", "fmn",
    "hydrogencarbonate", "pyrophosphate",
}


def _rhea_tsv(query: str, *, columns: str, limit: int = 200,
              timeout: int = 40) -> List[List[str]]:
    """Run a Rhea REST query and return the TSV rows (excluding the header)."""
    params = urllib.parse.urlencode(
        {"query": query, "columns": columns, "format": "tsv", "limit": limit})
    text = _http_get(f"{RHEA_REST}?{params}", timeout=timeout)
    rows = [ln.split("\t") for ln in text.splitlines() if ln.strip()]
    return rows[1:] if rows else []


def _parse_rhea_side(side: str) -> List[Tuple[float, str]]:
    """Parse one side of a Rhea equation into (coefficient, compound-name) pairs."""
    out: List[Tuple[float, str]] = []
    for term in side.split(" + "):
        term = term.strip()
        if not term:
            continue
        m = re.match(r"^(\d+)\s+(.*)$", term)      # a leading integer coefficient
        if m:
            coeff, name = float(m.group(1)), m.group(2).strip()
        else:
            coeff, name = 1.0, term
        # Drop the indefinite article Rhea uses for compound classes ("a ubiquinone").
        name = re.sub(r"^(a|an)\s+", "", name).strip()
        if name:
            out.append((coeff, name))
    return out


def build_rhea_pathway_db(
    target: str,
    *,
    expand_steps: int = 1,
    max_reactions: int = 300,
    time_budget: float = 45.0,
    max_frontier: int = 40,
) -> Tuple[cobra.Model, str, str]:
    """Build a focused universal model around ``target`` from Rhea (open licence).

    ``target`` is a compound name (e.g. "ethanol") or a ChEBI id (``CHEBI:16236``).
    Reactions mentioning the target are collected; ``expand_steps`` further rounds pull
    in reactions of the neighbouring (non-currency) compounds. Bounded by ``max_reactions``,
    a per-round frontier cap and a wall-clock ``time_budget`` so a hub compound can't hang
    the fetch. Metabolites are keyed by curated name so they reconcile with compounds
    already loaded from other databases.
    """
    target = target.strip()
    seed = target
    label = target
    cache_slug = re.sub(r"[^A-Za-z0-9]+", "_", target).strip("_")[:40] or "query"
    cache_name = f"rhea_{cache_slug}_s{expand_steps}_m{max_reactions}.json"
    cache_path = _cached_model_path(cache_name)
    cached = _load_built_model(cache_name)
    if cached is not None:
        return cached, f"Rhea: {label}", cache_path

    start = time.monotonic()
    deadline = start + time_budget
    columns = "rhea-id,equation,ec"
    collected: Dict[str, Tuple[str, str]] = {}     # rhea-id -> (equation, ec)
    frontier = {seed}
    seen: set = set()

    for _round in range(max(1, expand_steps + 1)):
        next_frontier: set = set()
        for term in sorted(frontier)[:max_frontier]:
            if term in seen or time.monotonic() > deadline:
                continue
            seen.add(term)
            try:
                rows = _rhea_tsv(term, columns=columns, limit=max_reactions)
            except DatabaseError:
                continue
            for row in rows:
                if len(collected) >= max_reactions:
                    break
                rid = row[0].strip()
                equation = row[1].strip() if len(row) > 1 else ""
                ec = row[2].strip() if len(row) > 2 else ""
                if not rid or "=" not in equation or rid in collected:
                    continue
                collected[rid] = (equation, ec)
                # Neighbours for the next round: non-currency participants.
                for side in equation.split("=", 1):
                    for _coeff, name in _parse_rhea_side(side):
                        if name.lower() not in _RHEA_CURRENCY:
                            next_frontier.add(name)
            if len(collected) >= max_reactions:
                break
        if len(collected) >= max_reactions or time.monotonic() > deadline:
            break
        frontier = next_frontier - seen

    if not collected:
        raise DatabaseError(
            f"Rhea has no reactions matching '{target}'. Try a ChEBI id (CHEBI:…) or a "
            "different spelling, or use another database source.")

    model = cobra.Model("rhea_pathway_db")
    mets: Dict[str, cobra.Metabolite] = {}

    def _met(name: str) -> cobra.Metabolite:
        key = name.lower()
        if key not in mets:
            slug = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")[:40] or f"m{len(mets)}"
            m = cobra.Metabolite(f"{slug}_c", name=name, compartment="c")
            mets[key] = m
        return mets[key]

    reactions = []
    for rid, (equation, ec) in collected.items():
        lhs, rhs = equation.split("=", 1)
        coeffs: Dict[cobra.Metabolite, float] = {}
        for coeff, name in _parse_rhea_side(lhs):
            met = _met(name)
            coeffs[met] = coeffs.get(met, 0.0) - abs(coeff)
        for coeff, name in _parse_rhea_side(rhs):
            met = _met(name)
            coeffs[met] = coeffs.get(met, 0.0) + abs(coeff)
        coeffs = {m: c for m, c in coeffs.items() if abs(c) > 1e-12}
        if not coeffs:
            continue
        safe_id = rid.replace(":", "_")
        r = cobra.Reaction(safe_id, name=rid)
        r.lower_bound, r.upper_bound = -1000.0, 1000.0   # Rhea reactions are directionless
        r.annotation["rhea"] = rid
        r.add_metabolites(coeffs)
        if ec:
            r.annotation["ec-code"] = ec.split(";")[0].strip()
        reactions.append(r)
    model.add_reactions(reactions)
    _save_built_model(model, cache_name)
    return model, f"Rhea: {label}", cache_path


# --------------------------------------------------------------------------- #
# ModelSEED biochemistry (downloaded universal)
# --------------------------------------------------------------------------- #
_SEED_COMPARTMENTS = {"0": "c", "1": "e", "2": "p"}


def _download_modelseed_tables(force: bool = False) -> Tuple[str, str]:
    """Download (and cache) the ModelSEED reactions + compounds TSVs."""
    rxn = os.path.join(_cache_dir(), "modelseed_reactions.tsv")
    cpd = os.path.join(_cache_dir(), "modelseed_compounds.tsv")
    if force or not os.path.exists(rxn) or os.path.getsize(rxn) < 1_000_000:
        _download(MODELSEED_REACTIONS_URL, rxn)
    if force or not os.path.exists(cpd) or os.path.getsize(cpd) < 1_000_000:
        _download(MODELSEED_COMPOUNDS_URL, cpd)
    return rxn, cpd


def _parse_seed_stoichiometry(field: str) -> Optional[Dict[str, float]]:
    """Parse a ModelSEED ``stoichiometry`` field into ``{met_id: coeff}``.

    Format: ``coeff:cpdid:compartment:community:"name"`` terms joined by ``;``.
    """
    if not field or field in ("null", "nan"):
        return None
    coeffs: Dict[str, float] = {}
    for term in field.split(";"):
        parts = term.split(":")
        if len(parts) < 3:
            continue
        try:
            coeff = float(parts[0])
        except ValueError:
            continue
        cpd, comp = parts[1], parts[2]
        if not cpd.startswith("cpd"):
            continue
        met_id = f"{cpd}_{_SEED_COMPARTMENTS.get(comp, 'c')}"
        coeffs[met_id] = coeffs.get(met_id, 0.0) + coeff
    return coeffs or None


def build_modelseed_universal(only_balanced: bool = True) -> cobra.Model:
    """Build a cobra universal reaction model from the ModelSEED biochemistry.

    Downloads the ModelSEED ``reactions.tsv`` + ``compounds.tsv`` (cached), keeps
    non-obsolete reactions (mass/charge balanced when ``only_balanced``), and
    tags metabolites/reactions with their SEED ids + EC numbers so they match a
    host model across namespaces. The built model is cached for instant reuse.
    """
    import pandas as pd

    cache_name = f"modelseed_universal{'_balanced' if only_balanced else ''}.json"
    cached = _load_built_model(cache_name)
    if cached is not None:
        return cached

    rxn_path, cpd_path = _download_modelseed_tables()
    cpd_df = pd.read_csv(cpd_path, sep="\t", dtype=str, low_memory=False).fillna("")
    cpd_info: Dict[str, dict] = {}
    for _, row in cpd_df.iterrows():
        cpd_info[row["id"]] = {
            "name": row.get("name", "") or row["id"],
            "formula": row.get("formula", "") or "",
            "charge": row.get("charge", "") or "",
        }

    rxn_df = pd.read_csv(rxn_path, sep="\t", dtype=str, low_memory=False).fillna("")
    model = cobra.Model("modelseed_universal")
    mets: Dict[str, cobra.Metabolite] = {}

    def _met(met_id: str) -> cobra.Metabolite:
        if met_id not in mets:
            base = met_id.rsplit("_", 1)[0]
            info = cpd_info.get(base, {})
            comp = met_id.rsplit("_", 1)[-1]
            m = cobra.Metabolite(met_id, name=info.get("name", base),
                                 formula=info.get("formula", "") or None, compartment=comp)
            try:
                if info.get("charge") not in ("", None):
                    m.charge = int(float(info["charge"]))
            except (ValueError, TypeError):
                pass
            m.annotation["seed.compound"] = base
            mets[met_id] = m
        return mets[met_id]

    reactions = []
    for _, row in rxn_df.iterrows():
        if str(row.get("is_obsolete", "0")) == "1":
            continue
        if only_balanced and row.get("status", "") != "OK":
            continue
        coeffs = _parse_seed_stoichiometry(row.get("stoichiometry", ""))
        if not coeffs:
            continue
        rid = row["id"]
        r = cobra.Reaction(rid, name=row.get("name", "") or rid)
        direction = row.get("reversibility", "") or row.get("direction", "=")
        if direction == ">":
            r.bounds = (0.0, 1000.0)
        elif direction == "<":
            r.bounds = (-1000.0, 0.0)
        else:
            r.bounds = (-1000.0, 1000.0)
        r.add_metabolites({_met(mid): c for mid, c in coeffs.items()})
        r.annotation["seed.reaction"] = rid
        ec = row.get("ec_numbers", "")
        if ec and ec not in ("null", "nan"):
            r.annotation["ec-code"] = [e for e in re.split(r"[|;]", ec) if e]
        reactions.append(r)

    model.add_reactions(reactions)
    _save_built_model(model, cache_name)
    return model


# --------------------------------------------------------------------------- #
# MetaNetX (downloaded universal)
# --------------------------------------------------------------------------- #
def download_metanetx_reactions(force: bool = False) -> str:
    """Download (and cache) MetaNetX ``reac_prop.tsv`` (~10 MB). Returns its path."""
    dest = os.path.join(_cache_dir(), "metanetx_reac_prop.tsv")
    if force or not os.path.exists(dest) or os.path.getsize(dest) < 1_000_000:
        _download(METANETX_REAC_URL, dest)
    return dest


def metanetx_chem_xref_path() -> str:
    return os.path.join(_cache_dir(), "metanetx_chem_xref.tsv")


def metanetx_reference_available() -> bool:
    """True if the (large) MetaNetX compound cross-reference file is downloaded."""
    p = metanetx_chem_xref_path()
    return os.path.exists(p) and os.path.getsize(p) > 100_000_000


def download_metanetx_chem_xref(force: bool = False) -> str:
    """Download (and cache) MetaNetX ``chem_xref.tsv`` (~680 MB, one-time).

    This file maps every external compound id (KEGG, BiGG, ChEBI, SEED, MetaCyc,
    HMDB, …) to an MNXM id and carries the compound name, so MetaNetX metabolites
    can be named and cross-referenced against any model. Stored permanently.
    """
    dest = metanetx_chem_xref_path()
    if force or not metanetx_reference_available():
        _download(METANETX_CHEM_XREF_URL, dest)
    return dest


def _parse_xref_source(source: str) -> Optional[Tuple[str, str]]:
    """Map a chem_xref ``source`` (``prefix:value``) to (annotation_key, value)."""
    if ":" not in source:
        return None
    prefix, value = source.split(":", 1)
    ns = _XREF_SOURCE_NS.get(prefix.lower())
    if not ns or not value:
        return None
    value = value.strip()
    if value.startswith("M_"):          # model-form ids: M_oh1 -> oh1
        value = value[2:]
    if ns == "chebi" and not value.upper().startswith("CHEBI:"):
        value = "CHEBI:" + value
    return ns, value


def build_metanetx_chem_index(used_ids: set, path: Optional[str] = None) -> Dict[str, dict]:
    """Stream chem_xref once, returning ``{MNXM: {"name", "ann": {key: [vals]}}}``
    restricted to the metabolite ids actually used (keeps it small)."""
    path = path or metanetx_chem_xref_path()
    index: Dict[str, dict] = {}
    if not os.path.exists(path):
        return index
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            source, mnx = parts[0], parts[1]
            if mnx not in used_ids:
                continue
            entry = index.setdefault(mnx, {"name": "", "ann": {}})
            desc = parts[2] if len(parts) > 2 else ""
            if desc and not entry["name"]:
                entry["name"] = desc.split("||")[0].strip()
            parsed = _parse_xref_source(source)
            if parsed:
                key, value = parsed
                vals = entry["ann"].setdefault(key, [])
                if value not in vals:
                    vals.append(value)
    return index


_MNX_TERM_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s+(\S+?)(?:@(\S+))?\s*$")


def _parse_mnx_side(side: str) -> List[Tuple[float, str, str]]:
    terms = []
    for raw in side.split(" + "):
        raw = raw.strip()
        if not raw:
            continue
        m = _MNX_TERM_RE.match(raw)
        if not m:
            # term without explicit coefficient
            parts = raw.split("@")
            terms.append((1.0, parts[0], parts[1] if len(parts) > 1 else "MNXD1"))
            continue
        coef, met, comp = m.group(1), m.group(2), m.group(3) or "MNXD1"
        terms.append((float(coef), met, comp))
    return terms


def build_metanetx_universal(
    *,
    only_balanced: bool = True,
    only_enzymatic: bool = False,
    max_reactions: Optional[int] = None,
    force_download: bool = False,
    path: Optional[str] = None,
    enrich: bool = False,
    xref_path: Optional[str] = None,
) -> cobra.Model:
    """Build a COBRA universal model from the MetaNetX ``reac_prop`` table.

    ``only_balanced`` keeps only mass/charge-balanced reactions (recommended —
    drops generic/polymer templates). ``only_enzymatic`` further restricts to
    reactions carrying an EC number. Transport reactions are always skipped (a
    universal pathway database describes chemistry, not membrane transport).

    With ``enrich`` (and the chem_xref reference file available/downloadable),
    metabolites get their real names and cross-reference annotations (KEGG, BiGG,
    ChEBI, …), so they can be matched against any model and read in plain language.
    """
    # Reuse a previously built+cached model unless the caller forces a refresh
    # or supplies a custom file (tests) or a reaction cap.
    cache_name = ""
    if path is None and max_reactions is None:
        kind = "enzymatic" if only_enzymatic else ("balanced" if only_balanced else "all")
        # '_named2' includes reaction names/short-ids (reac_xref); the bump forces a
        # one-time rebuild for users whose cached '_named' model predates that.
        cache_name = f"metanetx_universal_{kind}{'_named2' if enrich else ''}.json"
        if not force_download:
            cached = _load_built_model(cache_name)
            if cached is not None:
                return cached

    if path is None:
        path = download_metanetx_reactions(force=force_download)
    model = cobra.Model("metanetx_universal")
    mets: Dict[str, cobra.Metabolite] = {}
    used_mnx: set = set()

    def _met(mnx: str, comp: str) -> cobra.Metabolite:
        # Generic compartment MNXD1 -> "c"; keep others as a short suffix.
        c = "c" if comp in ("MNXD1", "") else comp.replace("MNXD", "d")
        mid = f"{mnx}_{c}"
        if mid not in mets:
            mets[mid] = cobra.Metabolite(mid, name=mnx, compartment=c)
            used_mnx.add(mnx)
        return mets[mid]

    reactions: List[cobra.Reaction] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            rid, equation, _ref, classifs, balanced, transport = cols[:6]
            if rid == "EMPTY" or transport.strip() == "T":
                continue
            if only_balanced and balanced.strip() != "B":
                continue
            ec = [e for e in classifs.split(";") if _EC_RE.fullmatch(e.strip())]
            if only_enzymatic and not ec:
                continue
            if " = " not in equation:
                continue
            lhs, _, rhs = equation.partition(" = ")
            try:
                left = _parse_mnx_side(lhs)
                right = _parse_mnx_side(rhs)
            except Exception:  # noqa: BLE001
                continue
            if not left and not right:
                continue
            r = cobra.Reaction(rid, name=rid)
            r.lower_bound, r.upper_bound = -1000.0, 1000.0
            coeffs: Dict[cobra.Metabolite, float] = {}
            for coef, mnx, comp in left:
                m = _met(mnx, comp)
                coeffs[m] = coeffs.get(m, 0.0) - coef
            for coef, mnx, comp in right:
                m = _met(mnx, comp)
                coeffs[m] = coeffs.get(m, 0.0) + coef
            coeffs = {m: c for m, c in coeffs.items() if c != 0.0}
            if not coeffs:
                continue
            r.add_metabolites(coeffs)
            if ec:
                r.annotation["ec-code"] = ec
            reactions.append(r)
            if max_reactions and len(reactions) >= max_reactions:
                break
    model.add_reactions(reactions)

    if enrich:
        if xref_path is None:
            xref_path = download_metanetx_chem_xref()
        _enrich_metanetx_metabolites(model, used_mnx, xref_path)
        try:
            used_rids = {r.id for r in model.reactions}
            _enrich_metanetx_reactions(model, used_rids, download_metanetx_reac_xref())
        except Exception:  # noqa: BLE001 - reaction names are a nicety, never fatal
            pass

    if cache_name:
        _save_built_model(model, cache_name)
    return model


def _enrich_metanetx_metabolites(model: cobra.Model, used_mnx: set, xref_path: str) -> None:
    """Attach names and cross-reference annotations from chem_xref to MNXM metabolites."""
    index = build_metanetx_chem_index(used_mnx, xref_path)
    for met in model.metabolites:
        base = met.id.rsplit("_", 1)[0]
        entry = index.get(base)
        if not entry:
            continue
        if entry.get("name"):
            met.name = entry["name"]
        ann = dict(met.annotation or {})
        ann.setdefault("metanetx.chemical", base)
        for key, vals in entry["ann"].items():
            ann[key] = vals
        met.annotation = ann


def metanetx_reac_xref_path() -> str:
    return os.path.join(_cache_dir(), "metanetx_reac_xref.tsv")


def download_metanetx_reac_xref(force: bool = False) -> str:
    """Download (and cache) MetaNetX ``reac_xref.tsv`` (~80 MB): reaction names and
    external short ids (BiGG/KEGG/SEED/MetaCyc) for MNXR ids."""
    dest = metanetx_reac_xref_path()
    if force or not (os.path.exists(dest) and os.path.getsize(dest) > 10_000_000):
        _download(METANETX_REAC_XREF_URL, dest)
    return dest


def _looks_like_name(text: str, ext_id: str) -> bool:
    if not text or text == ext_id:
        return False
    return not any(sym in text for sym in ("=", "<?>", "@", "||"))


def build_metanetx_reac_index(used_rids: set, path: Optional[str] = None) -> Dict[str, dict]:
    """Stream reac_xref once -> ``{MNXR: {"name", "short", "ann": {key: [vals]}}}``
    for the reactions actually used."""
    path = path or metanetx_reac_xref_path()
    index: Dict[str, dict] = {}
    if not os.path.exists(path):
        return index
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            source, mnxr = parts[0], parts[1]
            if mnxr not in used_rids or ":" not in source:
                continue
            prefix, value = source.split(":", 1)
            ns = _REAC_XREF_SOURCE_NS.get(prefix.lower())
            if not ns or not value:
                continue
            entry = index.setdefault(mnxr, {"name": "", "short": "", "names": {}, "ann": {}})
            entry["ann"].setdefault(ns, [])
            if value not in entry["ann"][ns]:
                entry["ann"][ns].append(value)
            desc = parts[2] if len(parts) > 2 else ""
            first = desc.split("||")[0].strip() if desc else ""
            if _looks_like_name(first, value):
                entry["names"].setdefault(ns, first)
    # resolve a preferred name + short id per reaction
    for entry in index.values():
        for ns in _REAC_NAME_PRIORITY:
            if entry["names"].get(ns):
                entry["name"] = entry["names"][ns]
                break
        for ns in _REAC_SHORT_PRIORITY:
            if entry["ann"].get(ns):
                entry["short"] = entry["ann"][ns][0]
                break
        entry.pop("names", None)
    return index


def _enrich_metanetx_reactions(model: cobra.Model, used_rids: set, path: str) -> None:
    """Attach readable names + external ids (and a short id) to MNXR reactions."""
    index = build_metanetx_reac_index(used_rids, path)
    for rxn in model.reactions:
        entry = index.get(rxn.id)
        if not entry:
            continue
        if entry.get("name"):
            rxn.name = entry["name"]
        ann = dict(rxn.annotation or {})
        for key, vals in entry["ann"].items():
            ann[key] = vals
        if entry.get("short"):
            ann.setdefault("short_id", entry["short"])
        rxn.annotation = ann


# --------------------------------------------------------------------------- #
# UniProt (live)
# --------------------------------------------------------------------------- #
def uniprot_enzymes_for_ec(ec: str, *, reviewed_only: bool = True,
                           limit: int = 25) -> pd.DataFrame:
    """Return candidate enzymes for an EC number from UniProt.

    Columns: Accession, Protein, Organism, Length, Reviewed, URL.
    """
    ec = ec.strip()
    if not _EC_RE.fullmatch(ec) and not re.fullmatch(r"\d+(\.\d+){1,3}", ec):
        raise DatabaseError(f"'{ec}' is not a valid EC number (e.g. 4.2.3.27).")
    query = f"(ec:{ec})"
    if reviewed_only:
        query += " AND (reviewed:true)"
    params = urllib.parse.urlencode({
        "query": query,
        "fields": "accession,protein_name,organism_name,length,reviewed",
        "format": "tsv",
        "size": str(max(1, min(limit, 100))),
    })
    text = _http_get(f"{UNIPROT_REST}?{params}")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) <= 1:
        return pd.DataFrame(columns=["Accession", "Protein", "Organism", "Length",
                                     "Reviewed", "URL"])
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        while len(parts) < 5:
            parts.append("")
        acc, protein, organism, length, reviewed = parts[:5]
        rows.append({
            "Accession": acc,
            "Protein": protein,
            "Organism": organism,
            "Length": length,
            "Reviewed": reviewed,
            "URL": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
        })
    return pd.DataFrame(rows)
