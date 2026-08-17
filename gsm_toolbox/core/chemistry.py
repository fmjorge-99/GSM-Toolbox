"""Structure resolution and chemical sanity checks for designed reactions.

Two capabilities the rest of the toolbox builds on:

**1. Getting a structure for a metabolite that has none.** Genome-scale databases usually
store only identifiers, which blocks every structure-based feature (isomer checking,
reaction-similarity enzyme search, thermodynamics of novel compounds). ``metabolite_smiles``
first reads the annotation, and — when asked — falls back to fetching the structure online
from the metabolite's own cross-references (PubChem by InChIKey, KEGG by compound id),
caching the result on disk so a model is only ever resolved once.

**2. Catching chemically impossible reactions that pass every other check.** A reaction
converting **2-methylpropanal** (isobutyraldehyde) to **n-butanol** is mass- and
charge-balanced — both are C4H10O — so the balance checker accepts it and FBA reports a
healthy flux. It is nonetheless impossible: reducing a *branched* aldehyde cannot give a
*straight-chain* alcohol. ``skeletal_isomer_warnings`` compares the carbon skeletons of a
reaction's participants and flags exactly this class of database error.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import cache

# --------------------------------------------------------------------------------------
# 1. Structures
# --------------------------------------------------------------------------------------
_SMILES_CACHE: Optional[Dict[str, str]] = None
_MEM: Dict[str, str] = {}


def _cache_path() -> str:
    return os.path.join(cache.base_dir(), "structures", "smiles_cache.json")


#: Bumped when a resolution bug means the cached *negatives* are wrong. Successful
#: lookups stay valid — only the "nothing found" entries are discarded, because a
#: negative recorded by broken code would otherwise be permanent.
_CACHE_VERSION = 2
_VERSION_KEY = "__version__"


def _load_cache() -> Dict[str, str]:
    global _SMILES_CACHE
    if _SMILES_CACHE is None:
        try:
            with open(_cache_path()) as fh:
                _SMILES_CACHE = json.load(fh)
        except Exception:  # noqa: BLE001
            _SMILES_CACHE = {}
        if _SMILES_CACHE.get(_VERSION_KEY) != _CACHE_VERSION:
            # Identifier extraction used to miss "InChI Key"/"KEGG Compound" annotations
            # and identifiers.org URLs, so it fell back to a name search that included the
            # formula ("Butanal C4H8O") and cached the failure forever. Drop those.
            _SMILES_CACHE = {k: v for k, v in _SMILES_CACHE.items()
                             if v and k != _VERSION_KEY}
            _SMILES_CACHE[_VERSION_KEY] = _CACHE_VERSION
            _save_cache()
    return _SMILES_CACHE


def _save_cache() -> None:
    try:
        p = _cache_path()
        os.makedirs(os.path.dirname(p), exist_ok=True)
        tmp = p + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(_load_cache(), fh)
        os.replace(tmp, p)
    except Exception:  # noqa: BLE001
        pass


def annotated_smiles(met) -> str:
    """A SMILES taken purely from the metabolite's own annotation (no network)."""
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
    inchi = first("inchi", "InChI")
    if inchi:
        try:
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            m = Chem.MolFromInchi(inchi)
            if m is not None:
                return Chem.MolToSmiles(m)
        except Exception:  # noqa: BLE001
            pass
    return ""


def metabolite_smiles(met, *, online: bool = False) -> str:
    """SMILES for a cobra metabolite.

    Reads the annotation first. With ``online=True`` and no annotated structure, resolves
    it from the metabolite's cross-references (PubChem by InChIKey, then KEGG by compound
    id, then PubChem by name) and caches the answer, so the expensive lookup happens once
    per compound ever. Returns "" when nothing can be resolved.
    """
    smi = annotated_smiles(met)
    if smi or not online:
        return smi

    key = _cache_key(met)
    if not key:
        return ""
    if key in _MEM:
        return _MEM[key]
    store = _load_cache()
    if key in store:
        _MEM[key] = store[key]
        return store[key]

    found = _fetch_smiles(met)
    _MEM[key] = found
    store[key] = found          # cache negatives too: do not retry a hopeless lookup
    _save_cache()
    return found


def _cache_key(met) -> str:
    from ..gui.widgets import structure_fetcher as sf  # local: avoids a GUI import cycle
    try:
        _s, _i, ik, kegg, chebi = sf.metabolite_structure_hints(met)
    except Exception:  # noqa: BLE001
        return ""
    for prefix, val in (("ik", ik), ("kegg", kegg), ("chebi", chebi)):
        if val:
            return f"{prefix}:{val}"
    name = _clean_name(met).lower()
    return f"name:{name}" if name else ""


def _fetch_smiles(met) -> str:
    """Resolve a structure online from a metabolite's cross-references."""
    from ..gui.widgets import structure_fetcher as sf
    try:
        _s, _i, ik, kegg, _chebi = sf.metabolite_structure_hints(met)
    except Exception:  # noqa: BLE001
        return ""
    # InChIKey is unambiguous, so try it first.
    if ik:
        try:
            smi = sf._pubchem_smiles(inchikey=ik)
            if smi:
                return smi
        except Exception:  # noqa: BLE001
            pass
    if kegg:
        try:
            mol = sf._kegg_molblock(kegg)
            if mol:
                from rdkit import Chem, RDLogger
                RDLogger.DisableLog("rdApp.*")
                m = Chem.MolFromMolBlock(mol)
                if m is not None:
                    return Chem.MolToSmiles(m)
        except Exception:  # noqa: BLE001
            pass
    # Names are last: they collide (see the DMAP incident), so only use an unambiguous one.
    name = _clean_name(met)
    if name:
        try:
            if not sf._ambiguous_name(name):
                smi = sf._pubchem_smiles(name=name)
                if smi:
                    return smi
        except Exception:  # noqa: BLE001
            pass
    return ""


def _clean_name(met) -> str:
    """A metabolite's name with any trailing molecular formula removed.

    BiGG names routinely carry the formula as a suffix — *"Butanal C4H8O"*, *"ATP
    C10H12N5O13P3"*. Passed to a compound-name lookup that matches nothing, which is
    then cached as a permanent negative for a compound PubChem knows perfectly well.
    """
    name = (getattr(met, "name", "") or "").strip()
    if not name:
        return ""
    # Only strip a trailing token that is a plausible formula AND leaves a name behind.
    cleaned = re.sub(r"\s+([A-Z][a-z]?\d*){2,}$", "", name).strip()
    return cleaned or name


# --------------------------------------------------------------------------------------
# 2. Skeletal-isomer checking
# --------------------------------------------------------------------------------------
@dataclass
class IsomerFinding:
    """One suspicious substrate→product pair inside a reaction."""

    reaction_id: str
    substrate_id: str
    substrate_name: str
    product_id: str
    product_name: str
    formula: str
    detail: str = ""
    #: True when the skeletons were compared by name rather than by structure. Such a
    #: finding is a suspicion to check, not a proof — it must never silently drop a route.
    from_names: bool = False

    def as_sentence(self) -> str:
        prefix = "Possible" if self.from_names else "Impossible"
        return (f"{prefix} — {self.reaction_id}: "
                f"{self.substrate_name or self.substrate_id} and "
                f"{self.product_name or self.product_id} are both {self.formula} but have "
                f"different carbon skeletons — {self.detail}")


@dataclass
class IsomerReport:
    checked: int = 0                     # reactions we could actually check
    unresolved: int = 0                  # reactions lacking structures on both sides
    findings: List[IsomerFinding] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.findings

    @property
    def coverage(self) -> str:
        total = self.checked + self.unresolved
        return f"{self.checked}/{total}" if total else "0/0"


def carbon_framework(smiles: str) -> str:
    """The bare carbon skeleton of a molecule, as canonical SMILES.

    Every heteroatom is deleted and all bonds are set to single, leaving only how the
    carbons are connected to one another. This is the invariant that a reaction which
    *adds or removes atoms* must preserve:

    * butanal ``CCCC=O`` → ``CCCC`` (straight C4)
    * 2-methylpropanal ``CC(C)C=O`` → ``CC(C)C`` (branched C4)

    Reducing an aldehyde cannot turn a branched skeleton into a straight one, so a
    reaction whose framework changes while atoms are added/removed is impossible.
    """
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return ""
        ed = Chem.RWMol(m)
        for atom in reversed(list(ed.GetAtoms())):      # reversed: indices stay valid
            if atom.GetSymbol() != "C":
                ed.RemoveAtom(atom.GetIdx())
        skel = ed.GetMol()
        if skel.GetNumAtoms() == 0:
            return ""
        for bond in skel.GetBonds():                    # ignore saturation differences
            bond.SetBondType(Chem.BondType.SINGLE)
        for atom in skel.GetAtoms():
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
            atom.SetFormalCharge(0)
            atom.SetIsAromatic(False)
        Chem.SanitizeMol(skel, catchErrors=True)
        return Chem.MolToSmiles(skel)
    except Exception:  # noqa: BLE001
        return ""


def _skeleton(smiles: str) -> str:
    """InChIKey connectivity block — full-molecule identity, used for isomerisations."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return ""
        return Chem.MolToInchiKey(m).split("-", 1)[0]
    except Exception:  # noqa: BLE001
        return ""


def _branching_signature(smiles: str) -> Optional[str]:
    """A crude description of the carbon backbone: 'branched' or 'straight'."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        for atom in m.GetAtoms():
            if atom.GetSymbol() != "C":
                continue
            heavy_c = sum(1 for nb in atom.GetNeighbors() if nb.GetSymbol() == "C")
            if heavy_c >= 3:
                return "branched"
        return "straight"
    except Exception:  # noqa: BLE001
        return None


#: Nomenclature that states the carbon skeleton outright. Chemical names are systematic
#: about this: a compound is called ``n-`` because it is straight, and ``iso``/``neo``/
#: ``2-methyl`` because it is not. Where no structure can be resolved, the name is the
#: only skeleton information available — and it is often enough.
_BRANCHED_NAME_HINTS = ("iso", "neo", "tert-", "tert ", "t-butyl", "sec-",
                        "2-methyl", "3-methyl", "2,2-dimethyl", "2,3-dimethyl",
                        "branched")
_STRAIGHT_NAME_HINTS = ("n-", "straight", "normal-", "unbranched")


def branching_from_name(name: str) -> Optional[str]:
    """'branched' / 'straight' inferred from a chemical name, or None if it says nothing.

    A fallback for the very common case of a database that carries identifiers but no
    structures. It is a heuristic and is reported as a suspicion, never as proof: a name
    can be wrong or idiosyncratic in a way a structure cannot.
    """
    text = (name or "").strip().lower()
    if not text:
        return None
    # Check branched first: "isobutanol" contains no straight hint, but a name like
    # "n-... isomer" should not be read as straight.
    if any(hint in text for hint in _BRANCHED_NAME_HINTS):
        return "branched"
    if text.startswith(_STRAIGHT_NAME_HINTS) or " n-" in text:
        return "straight"
    return None


def _formula_of(met) -> str:
    return (getattr(met, "formula", "") or "").strip()


def _carbon_count(met) -> int:
    try:
        return int((met.elements or {}).get("C", 0))
    except Exception:  # noqa: BLE001
        return 0


@dataclass
class BackfillReport:
    """What ``resolve_missing_chemistry`` managed to recover."""

    considered: int = 0          # metabolites that were missing a structure or a formula
    structures_added: int = 0
    formulas_added: int = 0
    unresolved: List[str] = field(default_factory=list)   # names we could not resolve

    @property
    def anything_added(self) -> bool:
        return bool(self.structures_added or self.formulas_added)

    def sentence(self) -> str:
        if not self.considered:
            return ""
        bits = []
        if self.structures_added:
            bits.append(f"{self.structures_added} structure(s)")
        if self.formulas_added:
            bits.append(f"{self.formulas_added} formula(s)")
        if not bits:
            return (f"Looked up {self.considered} metabolite(s) with missing chemistry; "
                    "none could be resolved from their cross-references.")
        return (f"Automatically resolved {' and '.join(bits)} for {self.considered} "
                "metabolite(s) that had none.")


#: BiGG's own formula/charge table, harvested once by
#: ``tools/fetch_bigg_metabolite_chem.py``. Consulted before anything else because a
#: formula must come from the same source as the stoichiometry — a MetaNetX formula
#: against BiGG stoichiometry leaves ~1500 reactions differing by protons alone.
_BIGG_CHEM: Optional[Dict[str, dict]] = None


def _bigg_chem_table() -> Dict[str, dict]:
    global _BIGG_CHEM
    if _BIGG_CHEM is None:
        path = os.path.join(cache.databases_dir(), "bigg_metabolite_chem.json")
        try:
            with open(path, encoding="utf-8") as fh:
                loaded = json.load(fh)
            _BIGG_CHEM = loaded if isinstance(loaded, dict) else {}
        except Exception:  # noqa: BLE001 — the table is an optimisation, never required
            _BIGG_CHEM = {}
    return _BIGG_CHEM


def _bigg_base_id(met) -> str:
    """The compartment-free BiGG id (``atp_c`` → ``atp``), preferring the annotation."""
    ann = getattr(met, "annotation", None) or {}
    for key in ("bigg.metabolite", "bigg", "BiGG"):
        val = ann.get(key)
        if val:
            v = (val[0] if isinstance(val, (list, tuple)) else str(val)).strip()
            if v:
                return v.rsplit("_", 1)[0] if re.search(r"_[a-z]\d?$", v) else v
    mid = str(getattr(met, "id", "") or "")
    return re.sub(r"_[a-z]\d?$", "", mid)


def _bigg_api_chem(bigg_id: str) -> Tuple[str, Optional[int]]:
    """``(formula, charge)`` for one BiGG id from the BiGG API; ``("", None)`` if unknown.

    Only an unambiguous answer is accepted: several formulae means BiGG holds conflicting
    values across models, and guessing one into a balance check is worse than declining.
    """
    import urllib.request

    url = f"http://bigg.ucsd.edu/api/v2/universal/metabolites/{bigg_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "GSM_ToolBox"})
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read())
    except Exception:  # noqa: BLE001 — offline, 404, or a transient failure
        return "", None
    formulae = data.get("formulae") or []
    charges = data.get("charges") or []
    formula = str(formulae[0]).strip() if len(formulae) == 1 else ""
    charge = charges[0] if len(charges) == 1 else None
    return formula, charge


def fetch_formula(met, *, online: bool = True) -> Tuple[str, Optional[int], str]:
    """``(formula, charge, source)`` for a metabolite that has none.

    Tries, in order of how well the answer matches the stoichiometry it will be checked
    against: the harvested BiGG table, the BiGG API, then the metabolite's structure
    (annotated, or resolved from its InChIKey/KEGG/ChEBI cross-references) converted to a
    formula. ``source`` is "" when nothing could be resolved, and ``charge`` is ``None``
    whenever the source did not give an unambiguous one — an absent charge stays absent
    rather than becoming a fabricated 0.
    """
    base = _bigg_base_id(met)
    if base:
        entry = _bigg_chem_table().get(base) or {}
        formula = str(entry.get("formula") or "").strip()
        if formula:
            return formula, entry.get("charge"), "BiGG (cached)"
    if online and base:
        formula, charge = _bigg_api_chem(base)
        if formula:
            _bigg_chem_table()[base] = {"formula": formula, "charge": charge}
            return formula, charge, "BiGG API"
    smiles = metabolite_smiles(met, online=online)
    if smiles:
        formula = _formula_from_smiles(smiles)
        if formula:
            # A structure gives no charge; leaving it None keeps the charge check honest.
            return formula, None, "structure (PubChem/KEGG/ChEBI)"
    return "", None, ""


@dataclass
class FormulaFetchReport:
    """What :func:`fetch_missing_formulas` recovered."""

    considered: int = 0
    filled: List[str] = field(default_factory=list)      # "name (id) → C40H56O6"
    unresolved: List[str] = field(default_factory=list)  # "name (id)"
    sources: Dict[str, str] = field(default_factory=dict)  # metabolite id -> source

    @property
    def anything_added(self) -> bool:
        return bool(self.filled)

    def sentence(self) -> str:
        if not self.considered:
            return "Every participant already has a formula — nothing to fetch."
        if not self.filled:
            return (f"Could not resolve a formula for any of the {self.considered} "
                    "metabolite(s) that lack one: "
                    + ", ".join(self.unresolved[:4])
                    + (" …" if len(self.unresolved) > 4 else "") + ".")
        text = (f"Fetched {len(self.filled)} of {self.considered} missing formula(s): "
                + "; ".join(self.filled[:6]) + (" …" if len(self.filled) > 6 else "") + ".")
        if self.unresolved:
            text += (f" Still unresolved: {', '.join(self.unresolved[:4])}"
                     + (" …" if len(self.unresolved) > 4 else "") + ".")
        return text


def fetch_missing_formulas(model, metabolite_ids: List[str], *,
                           online: bool = True) -> FormulaFetchReport:
    """Fill in the formulas of the named metabolites, in place.

    This is the on-demand counterpart to the "cannot be checked — X has no formula"
    verdict: rather than telling the user which database gap defeated the check, close
    it. Metabolites that already have a formula are skipped, and one that cannot be
    resolved is reported as still missing — never given a made-up formula, which would
    turn an honest "unverifiable" into a confident wrong answer.
    """
    from . import balancing

    rep = FormulaFetchReport()
    seen = set()
    for mid in metabolite_ids:
        if mid in seen or not model.metabolites.has_id(mid):
            continue
        seen.add(mid)
        met = model.metabolites.get_by_id(mid)
        if _formula_of(met):
            continue
        rep.considered += 1
        label = balancing.metabolite_label(met)
        formula, charge, source = fetch_formula(met, online=online)
        if not formula:
            rep.unresolved.append(label)
            continue
        met.formula = formula
        if charge is not None and getattr(met, "charge", None) is None:
            met.charge = charge
        rep.filled.append(f"{label} → {formula}")
        rep.sources[mid] = source
    return rep


def missing_formula_metabolites(model, reaction_ids: List[str]) -> List[str]:
    """Ids of every participant of ``reaction_ids`` that has no usable formula."""
    from . import balancing

    out: List[str] = []
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            continue
        for met in balancing.missing_formula_metabolites(model.reactions.get_by_id(rid)):
            if met.id not in out:
                out.append(met.id)
    return out


def resolve_missing_chemistry(model, reaction_ids: List[str], *,
                              online: bool = True) -> BackfillReport:
    """Fill in the structures and formulas the checks need, in place.

    The isomer check needs a structure and the balance check needs a formula; a database
    that supplies neither makes both report "cannot be verified". Rather than telling the
    user to go and fetch the data, do it: resolve each metabolite's structure from its
    cross-references (cached, so it costs nothing the second time) and derive the formula
    from that structure where the database has none.

    Returns a :class:`BackfillReport` describing what was recovered, so the caller can say
    what happened instead of what the user ought to do.
    """
    rep = BackfillReport()
    seen = set()
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            continue
        for met in model.reactions.get_by_id(rid).metabolites:
            if met.id in seen:
                continue
            seen.add(met.id)
            has_structure = bool(annotated_smiles(met))
            has_formula = bool(_formula_of(met))
            if has_structure and has_formula:
                continue
            rep.considered += 1
            smiles = metabolite_smiles(met, online=online)
            if smiles and not has_structure:
                # Write it back so every later check — isomer, ΔG, enzyme search — sees it.
                try:
                    met.annotation["SMILES"] = smiles
                    rep.structures_added += 1
                except Exception:  # noqa: BLE001 - annotation may be read-only
                    pass
            if not has_formula:
                # BiGG's own table first (same source as the stoichiometry), structure
                # second — see `fetch_formula`.
                formula, charge, _source = fetch_formula(met, online=online)
                if formula:
                    met.formula = formula
                    if charge is not None and getattr(met, "charge", None) is None:
                        met.charge = charge
                    rep.formulas_added += 1
                else:
                    rep.unresolved.append(getattr(met, "name", "") or met.id)
            elif not smiles and not has_structure:
                rep.unresolved.append(getattr(met, "name", "") or met.id)
    return rep


def _formula_from_smiles(smiles: str) -> str:
    """Molecular formula of a SMILES string, in the Hill-system form cobra expects."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:  # noqa: BLE001 - RDKit absent or structure unparseable
        return ""
    # RDKit appends the net charge (e.g. "C6H11O9P-2"); cobra's parser reads that as
    # elements and would choke, and the charge is carried separately anyway.
    return re.sub(r"[+-]\d*$", "", formula)


def check_reaction_isomers(rxn, *, online: bool = False) -> Tuple[int, int, List[IsomerFinding]]:
    """Look for substrate/product pairs that share a formula but not a skeleton.

    Returns ``(checked, unresolved, findings)``. Only carbon-containing organic pairs are
    considered, and cofactor-sized molecules are ignored — the pattern this catches is a
    main-chain compound being silently swapped for its isomer.
    """
    subs, prods = [], []
    for met, coeff in rxn.metabolites.items():
        if _carbon_count(met) < 2:
            continue                          # not an organic backbone
        if _is_cofactor(met):
            continue                          # NAD(P)(H), CoA … carry the chemistry, not it
        (subs if coeff < 0 else prods).append(met)
    if not subs or not prods:
        return 0, 0, []

    findings: List[IsomerFinding] = []
    checked = unresolved = 0
    for s in subs:
        ns = _carbon_count(s)
        for p in prods:
            # Only compare compounds of the same carbon count: those are the ones that
            # should be the *same molecule before and after*, modulo functional group.
            if _carbon_count(p) != ns:
                continue
            s_smi = metabolite_smiles(s, online=online)
            p_smi = metabolite_smiles(p, online=online)
            fs = carbon_framework(s_smi) if s_smi else ""
            fp = carbon_framework(p_smi) if p_smi else ""
            if not fs or not fp:
                # No structure on one side or both. Before giving up, ask the names:
                # a straight-chain compound paired with a branched one of the same
                # carbon count is the signature of a mislabelled database entry, and it
                # is visible without any structure at all.
                finding = _name_based_finding(rxn, s, p, ns)
                if finding is not None:
                    findings.append(finding)
                else:
                    unresolved += 1
                continue
            checked += 1
            if fs == fp:
                continue                      # backbone preserved — chemistry is fine
            # The backbones differ. That is legitimate only for a true isomerisation,
            # where nothing is added or removed (identical molecular formulas). If atoms
            # changed as well, the skeleton cannot also have rearranged in one step.
            f_s, f_p = _formula_of(s), _formula_of(p)
            if f_s and f_p and f_s == f_p:
                continue                      # genuine isomerase/mutase — allowed
            bs, bp = _branching_signature(s_smi), _branching_signature(p_smi)
            if bs and bp and bs != bp:
                detail = (f"{getattr(s, 'name', '') or s.id} is {bs} and "
                          f"{getattr(p, 'name', '') or p.id} is {bp}; a reaction that also "
                          "adds or removes atoms cannot rearrange the carbon backbone, so "
                          "this conversion is not possible as written")
            else:
                detail = ("their carbon backbones differ while atoms are also added or "
                          "removed, which no single step can do — most likely a "
                          "mislabelled compound")
            findings.append(IsomerFinding(
                reaction_id=rxn.id, substrate_id=s.id,
                substrate_name=getattr(s, "name", "") or s.id,
                product_id=p.id, product_name=getattr(p, "name", "") or p.id,
                formula=f"C{ns}", detail=detail))
    return checked, unresolved, findings


def _name_based_finding(rxn, s, p, carbons: int) -> Optional[IsomerFinding]:
    """A skeleton mismatch visible in the names alone, when structures are unavailable.

    Written after a route was accepted that reduced 2-methylpropanal straight to
    n-butanol. That reaction is mass- and charge-balanced — C4H8O + 2[H] = C4H10O — so
    every numeric check passes, and neither compound carried a structure, so the
    structural check never ran. The names said "branched" and "straight" the whole time.
    """
    s_name = getattr(s, "name", "") or s.id
    p_name = getattr(p, "name", "") or p.id
    bs, bp = branching_from_name(s_name), branching_from_name(p_name)
    if not bs or not bp or bs == bp:
        return None
    # A true isomerase legitimately rearranges the skeleton — but only if nothing else
    # changes. Identical formulas mean this could be one, so say nothing.
    f_s, f_p = _formula_of(s), _formula_of(p)
    if f_s and f_p and f_s == f_p:
        return None
    return IsomerFinding(
        reaction_id=rxn.id, substrate_id=s.id, substrate_name=s_name,
        product_id=p.id, product_name=p_name, formula=f"C{carbons}",
        detail=(f"the names say {s_name} is {bs} and {p_name} is {bp}; a reaction that "
                f"also adds or removes atoms cannot rearrange the carbon backbone, so "
                f"this conversion is not possible as written. Neither compound carries a "
                f"structure, so this is inferred from nomenclature — confirm before "
                f"discarding the route"),
        from_names=True)


_COFACTOR_HINTS = ("coa", "nad", "nadp", "fad", "fmn", "adenosyl", "atp", "adp", "amp",
                   "acyl-carrier", "acp", "thiamine", "biotin", "folate", "quinone",
                   "quinol", "cytochrome", "ferredoxin", "glutathione", "ubiquin")


def _is_cofactor(met) -> bool:
    name = (getattr(met, "name", "") or "").lower()
    mid = (getattr(met, "id", "") or "").lower()
    return any(h in name or h in mid for h in _COFACTOR_HINTS)


def check_pathway_isomers(model, reaction_ids: List[str], *, online: bool = False
                          ) -> IsomerReport:
    """Run the skeletal-isomer check over every reaction of a designed route."""
    rep = IsomerReport()
    for rid in reaction_ids:
        if not model.reactions.has_id(rid):
            continue
        try:
            checked, unresolved, found = check_reaction_isomers(
                model.reactions.get_by_id(rid), online=online)
        except Exception:  # noqa: BLE001 — never let a check break a search
            continue
        if found:
            rep.findings.extend(found)
            rep.checked += 1
        elif checked:
            rep.checked += 1
        elif unresolved:
            rep.unresolved += 1
    return rep
