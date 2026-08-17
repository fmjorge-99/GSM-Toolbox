"""Fetch/render a 2-D molecule structure image, using the persistent cache first.

Preferred output is a clean, monochrome (black-and-white) skeletal structure drawn
with **RDKit** from the compound's SMILES/InChI — crisp and free of coloured atom
dots. The SMILES is taken from the metabolite's annotation when present, otherwise
looked up once from PubChem. If RDKit or a structure string is unavailable we fall
back to PubChem's coloured 2-D PNG. Everything is cached on disk (keyed so B&W and
coloured variants don't clash) so repeat views are instant and offline.
"""

from __future__ import annotations

import re
import urllib.parse
import urllib.request

from PySide6.QtCore import QThread, Signal

from ...core import cache


def _http(url: str) -> bytes:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "GSM-ToolBox"})
        with urllib.request.urlopen(req, timeout=12) as resp:
            if resp.status == 200:
                return resp.read()
    except Exception:  # noqa: BLE001
        pass
    return b""


_NAME_CACHE: dict = {}


def _good_synonym(s: str) -> bool:
    """A human-readable common name — not a registry number, CID, or all-caps code."""
    s = (s or "").strip()
    if not s or len(s) > 60:
        return False
    low = s.lower()
    if low.startswith(("cid", "chebi", "chembl", "unii", "einecs", "ec ", "dsstox",
                       "schembl", "ac1", "mfcd", "ntp", "nsc", "hsdb", "brn ")):
        return False
    if any(ch.isdigit() for ch in s) and sum(ch.isdigit() for ch in s) > len(s) / 2:
        return False              # mostly digits = a registry id
    if "-" in s and s.replace("-", "").isdigit():
        return False              # CAS-like
    return any(ch.isalpha() for ch in s)


def name_from_inchikey(inchikey: str) -> str:
    """A common compound name for an InChIKey, from PubChem synonyms (cached).

    PubChem lists synonyms roughly by prevalence, so the first *readable* one is
    usually the everyday name. Falls back to the IUPAC name. "" if nothing resolves.
    """
    inchikey = (inchikey or "").strip()
    if not inchikey:
        return ""
    if inchikey in _NAME_CACHE:
        return _NAME_CACHE[inchikey]
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    name = ""
    data = _http(f"{base}/inchikey/{inchikey}/synonyms/JSON")
    if data:
        try:
            import json
            syns = [s for s in
                    json.loads(data)["InformationList"]["Information"][0]["Synonym"]
                    if _good_synonym(s)]
            # Prefer a name that already has lower-case letters (a proper common name)
            # over an ALL-CAPS registry-style synonym; title-case an all-caps fallback.
            name = next((s for s in syns if any(c.islower() for c in s)), "")
            if not name and syns:
                name = syns[0]
                if name.isupper():
                    name = name.capitalize()
        except Exception:  # noqa: BLE001
            name = ""
    if not name:
        d2 = _http(f"{base}/inchikey/{inchikey}/property/IUPACName/TXT")
        if d2:
            cand = d2.decode("utf-8", "ignore").strip().splitlines()
            if cand and cand[0] and len(cand[0]) <= 80:
                name = cand[0].strip()
    _NAME_CACHE[inchikey] = name
    return name


def name_from_smiles(smiles: str) -> str:
    """A common compound name for a SMILES, via its InChIKey (online). "" if none."""
    if not smiles:
        return ""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return ""
        return name_from_inchikey(Chem.MolToInchiKey(m))
    except Exception:  # noqa: BLE001
        return ""


def metabolite_structure_hints(met) -> tuple:
    """Return ``(smiles, inchi, inchikey, kegg, chebi)`` from a metabolite's
    annotation — the identifiers used to fetch a *correct* structure. KEGG/ChEBI
    ids are authoritative for metabolic compounds and avoid the name collisions
    that plague abbreviations (e.g. 'DMAP').

    Annotation keys are matched case- and punctuation-insensitively, and values are
    stripped of any identifiers.org URL. Both matter: BiGG-derived models write
    ``"InChI Key": "https://identifiers.org/inchikey/ZTQSA…"`` and ``"KEGG Compound"``,
    neither of which an exact-key lookup for ``inchi_key``/``kegg.compound`` finds. When
    that lookup silently returns nothing the caller falls back to searching by *name* —
    which for these models is something like ``"Butanal C4H8O"``, resolves to nothing,
    and gets cached as a permanent negative.
    """
    ann = getattr(met, "annotation", None) or {}
    # Normalise once: "InChI Key" / "inchi_key" / "InChIKey" all collapse to "inchikey".
    norm = {}
    for key, value in ann.items():
        slug = re.sub(r"[^a-z0-9]", "", str(key).lower())
        if slug not in norm:
            norm[slug] = value

    def _first(*slugs):
        for s in slugs:
            v = norm.get(s)
            if not v:
                continue
            raw = (v[0] if isinstance(v, (list, tuple)) else str(v)).strip()
            if not raw:
                continue
            if "://" in raw:            # identifiers.org URL → bare accession
                raw = raw.rstrip("/").rsplit("/", 1)[-1]
            return raw.strip()
        return ""

    chebi = _first("chebi", "chebiid")
    if chebi:
        chebi = chebi.upper().replace("CHEBI:", "").strip()
    kegg = _first("keggcompound", "kegg", "keggid", "keggligand")
    return (_first("smiles"),
            _first("inchi"),
            _first("inchikey"),
            kegg,
            chebi)


def _render_bw(smiles: str = "", inchi: str = "", molblock: str = "",
               size: int = 260) -> bytes:
    """Draw a clean black-and-white 2-D structure with RDKit; b"" if not possible."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception:  # noqa: BLE001 - RDKit missing -> caller falls back
        return b""
    mol = None
    try:
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
        if mol is None and inchi:
            mol = Chem.MolFromInchi(inchi, sanitize=True, removeHs=True)
        if mol is None and molblock:
            mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=True)
    except Exception:  # noqa: BLE001
        mol = None
    if mol is None:
        return b""
    try:
        d = rdMolDraw2D.MolDraw2DCairo(size, size)
        opts = d.drawOptions()
        opts.useBWAtomPalette()
        opts.clearBackground = True
        opts.bondLineWidth = 2
        opts.padding = 0.08
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:  # noqa: BLE001
        return b""


_COFACTOR_TAGS = [
    ("-coa", "CoA"), (" coa", "CoA"), ("-acp", "ACP"), ("-[acp]", "ACP"),
    ("-trna", "tRNA"), ("-glutathione", "GS"),
]


def split_cofactor(name: str) -> tuple:
    """If ``name`` is a core fused to a common carrier (acyl-CoA, -ACP…), return
    ``(core_name, tag)`` — e.g. ('malonyl-CoA' -> 'malonyl', 'CoA'); else (name, '')."""
    low = (name or "").lower()
    for suf, tag in _COFACTOR_TAGS:
        if low.endswith(suf):
            return name[: len(name) - len(suf)].strip(" -"), tag
    return name, ""


def _annotate_cofactor(png: bytes, tag: str, size: int) -> bytes:
    """Overlay a small '–<tag>' chip on a rendered core structure (e.g. '–CoA')."""
    try:
        from PySide6.QtCore import QRect, Qt
        from PySide6.QtGui import QColor, QFont, QImage, QPainter
    except Exception:  # noqa: BLE001
        return png
    img = QImage()
    if not img.loadFromData(png):
        return png
    p = QPainter(img)
    f = QFont(); f.setPointSize(max(9, size // 16)); f.setBold(True)
    p.setFont(f); p.setPen(QColor("#1a73e8"))
    p.drawText(QRect(4, img.height() - size // 6 - 4, img.width() - 8, size // 6),
               Qt.AlignRight | Qt.AlignVCenter, f"–{tag}")
    p.end()
    from PySide6.QtCore import QBuffer, QByteArray
    ba = QByteArray(); buf = QBuffer(ba); buf.open(QBuffer.WriteOnly)
    img.save(buf, "PNG")
    return bytes(ba)


def _pubchem_smiles(name: str = "", inchikey: str = "") -> str:
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    for path in ([f"inchikey/{inchikey}"] if inchikey else []) + \
                ([f"name/{urllib.parse.quote(name)}"] if name else []):
        data = _http(f"{base}/{path}/property/CanonicalSMILES/TXT")
        if data:
            s = data.decode("utf-8", "ignore").strip().splitlines()
            if s and s[0]:
                return s[0].strip()
    return ""


def _kegg_molblock(kegg_id: str) -> str:
    """Fetch the authoritative MOL block for a KEGG compound (e.g. C00235)."""
    kegg_id = (kegg_id or "").strip()
    if not kegg_id or kegg_id[0] not in "CGD":
        return ""
    data = _http(f"https://rest.kegg.jp/get/{kegg_id}/mol")
    return data.decode("utf-8", "ignore") if data else ""


def _ambiguous_name(name: str) -> bool:
    """True for short abbreviations / codes that resolve to the WRONG compound in a
    PubChem name search (e.g. 'DMAP' -> 4-dimethylaminopyridine, 'FPP', 'GPP').
    Such names must not be trusted for structure lookup."""
    n = (name or "").strip()
    if not n or " " in n:               # multi-word chemical names are usually specific
        return False
    compact = n.replace("-", "").replace(".", "")
    if len(compact) <= 6 and compact.isalnum():
        # mostly-uppercase short token = almost certainly an abbreviation
        letters = [c for c in compact if c.isalpha()]
        if letters and sum(c.isupper() for c in letters) >= max(2, len(letters) - 1):
            return True
    return False


def fetch_structure_png(name: str = "", inchikey: str = "", smiles: str = "",
                        inchi: str = "", kegg: str = "", chebi: str = "",
                        size: int = 220, bw: bool = True) -> bytes:
    """Synchronously return a *correct* 2-D structure PNG.

    Priority (most trustworthy first, avoiding name collisions like 'DMAP'):
    annotation SMILES/InChI → InChIKey → KEGG MOL → ChEBI → PubChem — and a plain
    name is used ONLY when it is not a short abbreviation. Cached on disk keyed by
    the most specific identifier so a wrong name lookup can't poison the cache.
    Safe off the GUI thread; never raises."""
    name = (name or "").strip()
    inchikey = (inchikey or "").strip()
    kegg = (kegg or "").strip()
    chebi = (chebi or "").strip()
    name_ok = bool(name) and not _ambiguous_name(name)

    # Cache key = the most specific identifier available (never an abbreviation).
    ident = (inchikey.upper() or (f"kegg:{kegg}" if kegg else "")
             or (f"chebi:{chebi}" if chebi else "") or (name.lower() if name_ok else ""))
    if not ident:
        ident = name.lower()          # last resort, may be an abbreviation
    key = f"struct::{'bw' if bw else 'color'}::{ident}::{size}"
    data = cache.get_image(key)
    if data:
        return data

    if bw:
        # 1) structure strings straight from the annotation
        data = _render_bw(smiles=smiles, inchi=inchi, size=max(size, 260))
        # 2) InChIKey -> PubChem canonical SMILES
        if not data and inchikey:
            data = _render_bw(smiles=_pubchem_smiles(inchikey=inchikey), size=max(size, 260))
        # 3) KEGG MOL (authoritative for metabolic compounds — fixes DMAP etc.)
        if not data and kegg:
            data = _render_bw(molblock=_kegg_molblock(kegg), size=max(size, 260))
        # 4) ChEBI id via PubChem xref
        if not data and chebi:
            s = _pubchem_smiles(name="") or ""
            xref = _http("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/"
                         f"RegistryID/CHEBI:{chebi}/property/CanonicalSMILES/TXT")
            if xref:
                lines = xref.decode("utf-8", "ignore").splitlines()
                if lines and lines[0].strip():
                    data = _render_bw(smiles=lines[0].strip(), size=max(size, 260))
        # 5) specific (non-abbreviation) name -> PubChem SMILES
        if not data and name_ok:
            data = _render_bw(smiles=_pubchem_smiles(name=name), size=max(size, 260))
    if not data and bw:
        # Cofactor conjugate (acyl-CoA, -ACP…): draw the core + a '–CoA' chip.
        core, tag = split_cofactor(name)
        if tag and core and not _ambiguous_name(core):
            core_png = _render_bw(smiles=_pubchem_smiles(name=core), size=max(size, 260))
            if core_png:
                data = _annotate_cofactor(core_png, tag, max(size, 260))
    if not data:                      # fall back to PubChem's coloured PNG (reliable ids only)
        base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
        px = f"?image_size={size}x{size}"
        if inchikey:
            data = _http(f"{base}/inchikey/{inchikey}/PNG{px}")
        if not data and name_ok:
            data = _http(f"{base}/name/{urllib.parse.quote(name)}/PNG{px}")
    if data:
        cache.put_image(key, data)
    return data


class StructureFetcher(QThread):
    fetched = Signal(str, bytes)  # tag (e.g. metabolite id), png bytes (empty on failure)

    def __init__(self, tag: str, name: str = "", inchikey: str = "", size: int = 220,
                 smiles: str = "", inchi: str = "", kegg: str = "", chebi: str = ""):
        super().__init__()
        self._tag = tag
        self._name = name or ""
        self._inchikey = inchikey or ""
        self._smiles = smiles or ""
        self._inchi = inchi or ""
        self._kegg = kegg or ""
        self._chebi = chebi or ""
        self._size = size

    @classmethod
    def for_metabolite(cls, met, size: int = 220):
        smiles, inchi, inchikey, kegg, chebi = metabolite_structure_hints(met)
        return cls(met.id, name=met.name or "", inchikey=inchikey, size=size,
                   smiles=smiles, inchi=inchi, kegg=kegg, chebi=chebi)

    def run(self) -> None:
        data = fetch_structure_png(name=self._name, inchikey=self._inchikey,
                                   smiles=self._smiles, inchi=self._inchi,
                                   kegg=self._kegg, chebi=self._chebi, size=self._size)
        self.fetched.emit(self._tag, data)
