"""Omics *dataset preparation*: turn a raw omics table into the two-column
``(model_id, value)`` format the toolbox's context-specific methods need, mapping
the data's identifiers onto the model's namespace. Pure-Python, no Qt.

Why this is needed — how omics data is actually distributed
----------------------------------------------------------
Public omics data almost never arrives as ``(model_gene_id, value)``:

* **Transcriptomics** (RNA-seq / microarray): a *matrix* of genes x samples with
  raw counts, TPM, FPKM or normalised intensities (GEO series-matrix files,
  featureCounts/HTSeq/Salmon/DESeq2 outputs, supplementary CSV/TSV/Excel). Gene
  ids are locus tags (``b0002``, ``sll1234``/``slr0009``), Ensembl/RefSeq ids or
  gene symbols. You usually must pick one sample/condition or average replicates.
* **Proteomics**: MaxQuant ``proteinGroups.txt`` (LFQ/iBAQ intensity columns),
  spectral counts, DIA reports. Ids are UniProt accessions (often several per row,
  ``P0A6;P0A7``, or FASTA headers ``sp|P0A6|…``), gene names, or locus tags.
* **Metabolomics**: a table of metabolites x samples with peak intensities. Ids are
  KEGG/HMDB/ChEBI ids, InChIKeys, or compound names — these map onto the model's
  *metabolites*, not its genes.

This module reads any of those (auto-detecting the delimiter / Excel), lets the
caller choose the id column, the value column(s) and how to aggregate, then maps
the ids onto the model via every identifier each gene/metabolite is known by
(id, name and cross-reference annotations). It returns the mapped ``{id: value}``
dict plus a :class:`MappingSummary` so the user can see how well it matched.
"""

from __future__ import annotations

import csv
import io
import re
import statistics
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import cobra
import pandas as pd

# omics kinds -> which model objects they map onto
GENE_KINDS = {"transcriptomics", "proteomics"}
METABOLITE_KINDS = {"metabolomics"}

# annotation keys that commonly hold a gene/protein identifier
_GENE_ANNOTATION_KEYS = (
    "ncbigene", "ncbiprotein", "refseq", "refseq_locus_tag", "uniprot",
    "kegg.genes", "kegg.gene", "ecogene", "asap", "locus_tag", "old_locus_tag",
    "protein_id", "gene", "genesymbol", "sgd", "ensembl",
)
_MET_ANNOTATION_KEYS = (
    "kegg.compound", "kegg.glycan", "chebi", "hmdb", "inchikey", "inchi_key",
    "bigg.metabolite", "metanetx.chemical", "seed.compound", "biocyc", "pubchem.compound",
)


class OmicsPrepError(Exception):
    """Raised when an omics table cannot be read or prepared."""


# --------------------------------------------------------------------------
# Reading a table robustly
# --------------------------------------------------------------------------
def read_table(path: str, *, max_preview_cols: int = 400) -> pd.DataFrame:
    """Read a CSV/TSV/Excel omics table, auto-detecting the delimiter. Comment
    lines (starting with ``#`` or ``!`` — the latter for GEO series-matrix files)
    are skipped. Returns a pandas DataFrame with string column labels."""
    low = path.lower()
    if low.endswith((".xlsx", ".xls")):
        try:
            df = pd.read_excel(path)
        except Exception as exc:  # noqa: BLE001
            raise OmicsPrepError(
                "Could not read the Excel file (is 'openpyxl' installed?):\n"
                f"{exc}") from exc
    else:
        # sniff delimiter from a sample of non-comment lines
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                lines = [ln for ln in fh.readlines()
                         if ln.strip() and not ln.lstrip().startswith(("#", "!"))]
        except Exception as exc:  # noqa: BLE001
            raise OmicsPrepError(f"Could not open the file:\n{exc}") from exc
        if not lines:
            raise OmicsPrepError("The file has no data rows.")
        sample = "".join(lines[:50])
        sep = "\t"
        try:
            sep = csv.Sniffer().sniff(sample, delimiters=",\t;| ").delimiter
        except Exception:  # noqa: BLE001
            # fall back to whichever common delimiter appears most in the header
            header = lines[0]
            sep = max(",\t;|", key=lambda d: header.count(d))
        try:
            df = pd.read_csv(io.StringIO("".join(lines)), sep=sep, engine="python")
        except Exception as exc:  # noqa: BLE001
            raise OmicsPrepError(f"Could not parse the table (delimiter '{sep}'):\n{exc}") from exc
    if df.shape[1] < 2:
        raise OmicsPrepError("Need at least two columns: an identifier column and a value column.")
    df.columns = [str(c).strip() for c in df.columns][:max_preview_cols] \
        if df.shape[1] <= max_preview_cols else [str(c).strip() for c in df.columns]
    return df


def numeric_columns(df: pd.DataFrame) -> List[str]:
    """Columns that look numeric (candidate value/sample columns)."""
    out = []
    for c in df.columns:
        col = pd.to_numeric(df[c], errors="coerce")
        if col.notna().mean() >= 0.6:      # mostly numbers
            out.append(str(c))
    return out


def guess_id_column(df: pd.DataFrame) -> str:
    """Best guess for the identifier column: a named id column, else the first
    non-numeric column, else the first column."""
    num = set(numeric_columns(df))
    named = re.compile(r"(gene|locus|protein|orf|id|name|symbol|accession|compound|"
                       r"metabolite|kegg|uniprot)", re.I)
    for c in df.columns:
        if str(c) not in num and named.search(str(c)):
            return str(c)
    for c in df.columns:
        if str(c) not in num:
            return str(c)
    return str(df.columns[0])


# --------------------------------------------------------------------------
# Identifier normalisation + model index
# --------------------------------------------------------------------------
def normalize_id(text) -> str:
    """Canonical form for matching: lower-cased, trimmed, without a trailing
    version (``.1``) or isoform (``-2``) suffix, quotes or surrounding brackets."""
    s = str(text).strip().strip('"\'')
    s = re.sub(r"\.\d+$", "", s)          # transcript/version suffix
    s = re.sub(r"-\d+$", "", s)           # protein isoform suffix
    return s.lower()


def _tokens(raw) -> List[str]:
    """Split a raw identifier cell into candidate ids: many omics tables pack
    several ids per cell (``P1;P2``), or a FASTA header (``sp|P12345|NAME_ORG``)."""
    s = str(raw).strip()
    if not s:
        return []
    # FASTA-style header: keep the middle accession
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", s, re.I)
    if m:
        return [m.group(1)]
    parts = re.split(r"[;,|\s]+", s)
    return [p for p in parts if p]


def _flatten_annotation(value) -> List[str]:
    if isinstance(value, (list, tuple, set)):
        out = []
        for v in value:
            out.extend(_flatten_annotation(v))
        return out
    return [str(value)]


def build_identifier_index(model: cobra.Model, kind: str) -> Dict[str, str]:
    """Map every identifier a gene (or metabolite) is known by -> its model id.

    For gene kinds this indexes each gene's id, name and gene-ish annotation
    values; for metabolomics it indexes metabolite id, name and chemical xrefs.
    Both raw and normalised forms are stored so matching is forgiving."""
    index: Dict[str, str] = {}

    def _add(key, target_id):
        if key is None:
            return
        for form in (str(key).strip(), normalize_id(key)):
            if form and form not in index:
                index[form] = target_id

    if kind in METABOLITE_KINDS:
        objs = model.metabolites
        ann_keys = _MET_ANNOTATION_KEYS
    else:
        objs = model.genes
        ann_keys = _GENE_ANNOTATION_KEYS

    for obj in objs:
        oid = obj.id
        _add(obj.id, oid)
        _add(getattr(obj, "name", None), oid)
        ann = getattr(obj, "annotation", None)
        if isinstance(ann, dict):
            for k in ann_keys:
                if k in ann:
                    for v in _flatten_annotation(ann[k]):
                        # chebi values may be "CHEBI:12345" -> also index bare number
                        _add(v, oid)
                        if ":" in v:
                            _add(v.split(":", 1)[1], oid)
    return index


def _lookup(raw, index: Dict[str, str]) -> Optional[str]:
    for tok in _tokens(raw):
        if tok in index:
            return index[tok]
        n = normalize_id(tok)
        if n in index:
            return index[n]
    return None


# --------------------------------------------------------------------------
# Aggregation + mapping
# --------------------------------------------------------------------------
_AGGREGATORS = {
    "mean": lambda xs: statistics.fmean(xs),
    "median": lambda xs: statistics.median(xs),
    "max": max,
    "min": min,
    "sum": sum,
    "first": lambda xs: xs[0],
}


@dataclass
class MappingSummary:
    kind: str
    id_column: str
    value_columns: List[str]
    aggregate: str
    n_source_rows: int = 0
    n_source_mapped: int = 0          # source rows that hit a model id
    n_model_targets: int = 0          # distinct model ids covered
    n_model_total: int = 0            # genes (or metabolites) in the model
    unmatched_examples: List[str] = field(default_factory=list)
    n_unmatched: int = 0

    @property
    def coverage(self) -> float:
        return (self.n_model_targets / self.n_model_total) if self.n_model_total else 0.0

    def text(self) -> str:
        tgt = "metabolites" if self.kind in METABOLITE_KINDS else "genes"
        lines = [
            f"Data type: {self.kind}",
            f"Identifier column: {self.id_column}",
            f"Value column(s): {', '.join(self.value_columns)}  (aggregate: {self.aggregate})",
            f"Source rows: {self.n_source_rows}  |  mapped to the model: {self.n_source_mapped}",
            f"Model {tgt} covered: {self.n_model_targets} / {self.n_model_total} "
            f"({self.coverage*100:.1f}%)",
            f"Unmatched identifiers: {self.n_unmatched}",
        ]
        if self.unmatched_examples:
            lines.append("Examples of unmatched ids: "
                         + ", ".join(self.unmatched_examples[:12]))
        return "\n".join(lines)


@dataclass
class PreparedDataset:
    kind: str
    values: Dict[str, float]          # {model_id: value}
    summary: MappingSummary
    unmatched: List[str] = field(default_factory=list)

    def to_frame(self) -> pd.DataFrame:
        col = "metabolite" if self.kind in METABOLITE_KINDS else "gene"
        return pd.DataFrame({col: list(self.values.keys()),
                             "value": list(self.values.values())})


def prepare_dataset(df: pd.DataFrame, model: cobra.Model, *, kind: str,
                    id_column: str, value_columns: List[str],
                    aggregate: str = "mean") -> PreparedDataset:
    """Map an omics table onto the model's namespace. ``kind`` is one of
    ``transcriptomics``/``proteomics`` (mapped to genes) or ``metabolomics``
    (mapped to metabolites). ``value_columns`` are aggregated per row (e.g. to
    average replicate samples), then re-aggregated across ids that collide."""
    if aggregate not in _AGGREGATORS:
        raise OmicsPrepError(f"Unknown aggregate '{aggregate}'.")
    if id_column not in df.columns:
        raise OmicsPrepError(f"No column '{id_column}' in the table.")
    value_columns = [c for c in value_columns if c in df.columns]
    if not value_columns:
        raise OmicsPrepError("Select at least one numeric value column.")

    index = build_identifier_index(model, kind)
    agg = _AGGREGATORS[aggregate]
    num = df[value_columns].apply(pd.to_numeric, errors="coerce")

    per_target: Dict[str, List[float]] = {}
    n_mapped = 0
    unmatched: List[str] = []
    ids = df[id_column].tolist()
    for i, raw in enumerate(ids):
        row_vals = [float(v) for v in num.iloc[i].tolist() if v == v]  # drop NaN
        if not row_vals:
            continue
        val = agg(row_vals)
        target = _lookup(raw, index)
        if target is None:
            if str(raw).strip():
                unmatched.append(str(raw).strip())
            continue
        per_target.setdefault(target, []).append(val)
        n_mapped += 1

    values = {tid: agg(vs) for tid, vs in per_target.items()}
    n_total = len(model.metabolites if kind in METABOLITE_KINDS else model.genes)
    # de-duplicate unmatched examples, keep order
    seen, uex = set(), []
    for u in unmatched:
        if u not in seen:
            seen.add(u)
            uex.append(u)
    summary = MappingSummary(
        kind=kind, id_column=id_column, value_columns=value_columns, aggregate=aggregate,
        n_source_rows=len(ids), n_source_mapped=n_mapped, n_model_targets=len(values),
        n_model_total=n_total, unmatched_examples=uex, n_unmatched=len(uex))
    return PreparedDataset(kind=kind, values=values, summary=summary, unmatched=uex)


def write_prepared_csv(dataset: PreparedDataset, path: str) -> None:
    """Save a prepared dataset as the two-column (id, value) CSV/TSV the toolbox
    loads for eFlux/GIMME."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    dataset.to_frame().to_csv(path, sep=sep, index=False)
