"""Gene and reaction deletion screens (essentiality analysis).

Wraps cobra's deletion functions and returns tidy, GUI-friendly DataFrames with
the growth ratio relative to wild-type, so essential genes/reactions are obvious.
"""

from __future__ import annotations

import os

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion

# On Windows, spawning worker processes (re-importing cobra + pickling the model
# per worker) costs far more than the per-deletion LP solve, so cobra's default
# multiprocessing is several times SLOWER than a single fast process. Force the
# single-process path on Windows; allow cobra's parallelism elsewhere (cheap fork).
_PROCESSES = 1 if os.name == "nt" else None


def _tidy(deletion_df: pd.DataFrame, wt_growth: float, kind: str) -> pd.DataFrame:
    """Flatten cobra's deletion result (frozenset index) into a flat table."""
    rows = []
    for _, row in deletion_df.iterrows():
        ids = row.get("ids")
        if isinstance(ids, (set, frozenset)):
            label = ", ".join(sorted(ids))
        else:
            label = str(ids)
        growth = row.get("growth", float("nan"))
        status = row.get("status", "")
        ratio = (growth / wt_growth) if wt_growth not in (0, None) else float("nan")
        # A knockout is essential if it abolishes growth: either the model becomes
        # infeasible, or growth drops below 1% of wild-type.
        is_nan = ratio != ratio
        essential = (status != "optimal") or (not is_nan and ratio < 0.01)
        rows.append({
            kind: label,
            "growth": growth,
            "growth_ratio": ratio,
            "status": status,
            "essential": bool(essential),
        })
    df = pd.DataFrame(rows).sort_values("growth_ratio", na_position="first")
    return df.reset_index(drop=True)


def single_reaction_deletions(model: cobra.Model) -> pd.DataFrame:
    """Knock out each reaction one at a time; report resulting growth."""
    wt = model.slim_optimize()
    result = single_reaction_deletion(model, processes=_PROCESSES)
    return _tidy(result, float(wt), "reaction")


def single_gene_deletions(model: cobra.Model) -> pd.DataFrame:
    """Knock out each gene one at a time; report resulting growth."""
    wt = model.slim_optimize()
    result = single_gene_deletion(model, processes=_PROCESSES)
    return _tidy(result, float(wt), "gene")
