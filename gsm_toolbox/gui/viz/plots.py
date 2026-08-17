"""Publication-grade plot builders for the graphical engine.

Every function takes a matplotlib ``Figure`` (cleared and drawn in place) plus its
data, and returns a tidy ``pandas.DataFrame`` of the values behind the figure — so
the shared canvas can offer "Export data (CSV)" alongside SVG/PDF/PNG for every
plot (the proposal's reproducibility requirement). Functions are defensive: on
empty/degenerate input they draw a short explanatory message instead of raising.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from . import theme


def _empty(fig, message: str) -> pd.DataFrame:
    fig.clear()
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True,
            color="#5f6368", fontsize=10)
    return pd.DataFrame()


# --- Strategy Explorer -------------------------------------------------------
def multi_strategy_heatmap(fig, matrix: Dict[str, List[float]], strategies: List[str],
                           labels: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Reactions × strategies flux heatmap (diverging), values annotated."""
    if not matrix or not strategies:
        return _empty(fig, "Add at least one strategy and some reactions to compare.")
    labels = labels or {}
    rxn_ids = list(matrix)
    data = np.array([matrix[r] for r in rxn_ids], dtype=float)
    fig.clear()
    ax = fig.add_subplot(111)
    vmax = np.nanmax(np.abs(data)) or 1.0
    im = ax.imshow(data, aspect="auto", cmap=theme.DIVERGING, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(strategies)))
    ax.set_xticklabels(strategies, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(rxn_ids)))
    ax.set_yticklabels([labels.get(r, r) for r in rxn_ids], fontsize=7)
    if data.shape[0] * data.shape[1] <= 200:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, f"{data[i, j]:.2g}", ha="center", va="center",
                        fontsize=6, color="#222")
    ax.set_title("Flux across strategies")
    fig.colorbar(im, ax=ax, label="flux (mmol gDW⁻¹ h⁻¹)")
    return pd.DataFrame(data, index=rxn_ids, columns=strategies).reset_index(names="reaction")


def titre_waterfall(fig, names: List[str], values: List[float], target: str = "") -> pd.DataFrame:
    """Product flux after each engineering round (a titre waterfall)."""
    if not names:
        return _empty(fig, "Save strategies with a product target to see the titre waterfall.")
    fig.clear()
    ax = fig.add_subplot(111)
    colors = [theme.UP if v >= 0 else theme.DOWN for v in values]
    ax.bar(range(len(names)), values, color=colors)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(f"|flux| of {target}" if target else "product flux")
    ax.set_title("Titre after each engineering round")
    theme.style_axes(ax)
    return pd.DataFrame({"strategy": names, "product_flux": values})


def parallel_coordinates(fig, matrix: Dict[str, List[float]], strategies: List[str],
                         labels: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """One line per reaction across the strategy axis (parallel coordinates)."""
    if not matrix or len(strategies) < 2:
        return _empty(fig, "Pick at least two strategies and some reactions.")
    labels = labels or {}
    fig.clear()
    ax = fig.add_subplot(111)
    x = range(len(strategies))
    for i, (rid, series) in enumerate(matrix.items()):
        ax.plot(x, series, marker="o", ms=3, lw=1.3,
                color=theme.CATEGORICAL[i % len(theme.CATEGORICAL)],
                label=labels.get(rid, rid))
    ax.set_xticks(list(x))
    ax.set_xticklabels(strategies, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("flux")
    ax.set_title("Flux trajectories across strategies")
    if len(matrix) <= 12:
        ax.legend(fontsize=7, ncol=2)
    theme.style_axes(ax)
    return pd.DataFrame(matrix, index=strategies).T.reset_index(names="reaction")


# --- Per-utility catalogue ---------------------------------------------------
def production_envelope(fig, df: pd.DataFrame, target: str,
                        overlays: Optional[Dict[str, pd.DataFrame]] = None) -> pd.DataFrame:
    """Filled feasible envelope; optionally overlay several strains/conditions."""
    if df is None or df.empty:
        return _empty(fig, "No production-envelope data.")
    fig.clear()
    ax = fig.add_subplot(111)
    x = df[target] if target in df.columns else df.iloc[:, 0]
    ax.fill_between(x, df["flux_minimum"], df["flux_maximum"], color=theme.CATEGORICAL[0],
                    alpha=0.22, label="feasible")
    ax.plot(x, df["flux_maximum"], color=theme.CATEGORICAL[0], lw=2)
    for i, (name, odf) in enumerate((overlays or {}).items(), start=1):
        ox = odf[target] if target in odf.columns else odf.iloc[:, 0]
        ax.plot(ox, odf["flux_maximum"], lw=1.8,
                color=theme.CATEGORICAL[i % len(theme.CATEGORICAL)], label=name)
    ax.set_xlabel(f"{target} flux")
    ax.set_ylabel("growth rate")
    ax.set_title("Production envelope")
    ax.legend(fontsize=8)
    theme.style_axes(ax)
    return df


def fva_tornado(fig, df: pd.DataFrame, top: int = 25) -> pd.DataFrame:
    """Whisker plot of [min, max] flux per reaction (widest ranges first)."""
    if df is None or df.empty or not {"minimum", "maximum"}.issubset(df.columns):
        return _empty(fig, "No FVA [min, max] data.")
    d = df.copy()
    rcol = "reaction" if "reaction" in d.columns else d.columns[0]
    d["_span"] = (d["maximum"] - d["minimum"]).abs()
    d = d.sort_values("_span", ascending=False).head(top).iloc[::-1]
    fig.clear()
    ax = fig.add_subplot(111)
    y = range(len(d))
    ax.hlines(y, d["minimum"], d["maximum"], color=theme.CATEGORICAL[0], lw=4, alpha=0.7)
    ax.plot(d["minimum"], y, "|", color=theme.DOWN, ms=8)
    ax.plot(d["maximum"], y, "|", color=theme.UP, ms=8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(d[rcol], fontsize=7)
    ax.axvline(0, color="#B8C0CC", lw=0.8)
    ax.set_xlabel("flux range")
    ax.set_title(f"FVA flux ranges (top {len(d)} by span)")
    theme.style_axes(ax)
    return df


def fseof_scan(fig, table: pd.DataFrame) -> pd.DataFrame:
    """Multi-line flux scan coloured by target class (amplify vs knock-down)."""
    needed = {"reaction", "flux_start", "flux_end", "target_type"}
    if table is None or table.empty or not needed.issubset(table.columns):
        return _empty(fig, "No FSEOF scan data.")
    fig.clear()
    ax = fig.add_subplot(111)
    for _, row in table.iterrows():
        up = "overexpression" in str(row["target_type"])
        ax.plot([0, 1], [row["flux_start"], row["flux_end"]],
                color=theme.UP if up else theme.DOWN, lw=1.3, alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["low enforced flux", "high enforced flux"], fontsize=8)
    ax.set_ylabel("reaction flux")
    ax.set_title("FSEOF scan (red = amplify, blue = knock-down)")
    theme.style_axes(ax)
    return table


def phase_plane(fig, df: pd.DataFrame, x_name: str = "x", y_name: str = "y") -> pd.DataFrame:
    """2-D heatmap of the objective over two exchange fluxes, with contours."""
    if df is None or df.empty or not {"x", "y", "objective"}.issubset(df.columns):
        return _empty(fig, "No phase-plane data.")
    pivot = df.pivot(index="y", columns="x", values="objective")
    fig.clear()
    ax = fig.add_subplot(111)
    im = ax.imshow(pivot.values, origin="lower", aspect="auto", cmap=theme.SEQUENTIAL,
                   extent=[pivot.columns.min(), pivot.columns.max(),
                           pivot.index.min(), pivot.index.max()])
    try:
        cs = ax.contour(pivot.columns, pivot.index, pivot.values, colors="white",
                        linewidths=0.6, alpha=0.7)
        ax.clabel(cs, fontsize=6, inline=True)
    except Exception:  # noqa: BLE001
        pass
    ax.set_xlabel(f"{x_name} flux")
    ax.set_ylabel(f"{y_name} flux")
    ax.set_title("Phenotype phase plane")
    fig.colorbar(im, ax=ax, label="objective")
    return df


def exchange_flux_bars(fig, fluxes: Dict[str, float], exchange_ids: List[str],
                       labels: Optional[Dict[str, str]] = None, top: int = 20) -> pd.DataFrame:
    """Uptake (negative) / secretion (positive) bar chart of exchange fluxes."""
    labels = labels or {}
    rows = [(rid, float(fluxes.get(rid, 0.0))) for rid in exchange_ids]
    rows = [r for r in rows if abs(r[1]) > 1e-9]
    if not rows:
        return _empty(fig, "No non-zero exchange fluxes in this state.")
    rows.sort(key=lambda r: abs(r[1]), reverse=True)
    rows = rows[:top][::-1]
    fig.clear()
    ax = fig.add_subplot(111)
    vals = [v for _, v in rows]
    colors = [theme.UP if v > 0 else theme.DOWN for v in vals]
    ax.barh(range(len(rows)), vals, color=colors)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([labels.get(r, r) for r, _ in rows], fontsize=7)
    ax.axvline(0, color="#B8C0CC", lw=0.8)
    ax.set_xlabel("flux  (– uptake / + secretion)")
    ax.set_title("Exchange fluxes")
    theme.style_axes(ax)
    return pd.DataFrame(rows, columns=["exchange", "flux"])


def strain_design_comparison(fig, df: pd.DataFrame) -> pd.DataFrame:
    """Grouped bars comparing designs: growth vs guaranteed vs max product."""
    cols = [c for c in ("predicted_growth", "guaranteed_product", "product_at_max_growth")
            if c in (df.columns if df is not None else [])]
    if df is None or df.empty or not cols:
        return _empty(fig, "No strain-design comparison data.")
    d = df.head(10)
    labels = d["knockouts"] if "knockouts" in d.columns else d.index.astype(str)
    fig.clear()
    ax = fig.add_subplot(111)
    n, w = len(cols), 0.8 / max(1, len(cols))
    x = np.arange(len(d))
    for i, c in enumerate(cols):
        ax.bar(x + i * w, pd.to_numeric(d[c], errors="coerce").fillna(0.0), w,
               label=c.replace("_", " "), color=theme.CATEGORICAL[i % len(theme.CATEGORICAL)])
    ax.set_xticks(x + w * (n - 1) / 2)
    ax.set_xticklabels([str(s)[:22] for s in labels], rotation=30, ha="right", fontsize=7)
    ax.set_title("Strain-design comparison")
    ax.legend(fontsize=8)
    theme.style_axes(ax)
    return df


def essentiality_bars(fig, df: pd.DataFrame, top: int = 25) -> pd.DataFrame:
    """Growth ratio per deletion (lowest = most essential)."""
    if df is None or df.empty:
        return _empty(fig, "No deletion data.")
    d = df.copy()
    gcol = next((c for c in ("growth", "objective_value", "growth_ratio") if c in d.columns), None)
    idcol = next((c for c in ("ids", "reaction", "gene", "id") if c in d.columns), d.columns[0])
    if gcol is None:
        return _empty(fig, "No growth column in the deletion result.")
    d = d.sort_values(gcol).head(top).iloc[::-1]
    fig.clear()
    ax = fig.add_subplot(111)
    ax.barh(range(len(d)), pd.to_numeric(d[gcol], errors="coerce").fillna(0.0),
            color=theme.SEQUENTIAL if False else theme.CATEGORICAL[2])
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([str(v)[:24] for v in d[idcol]], fontsize=7)
    ax.set_xlabel("mutant growth")
    ax.set_title(f"Deletion impact (lowest growth = most essential, top {len(d)})")
    theme.style_axes(ax)
    return df


def flux_sampling_violin(fig, samples: pd.DataFrame, top: int = 15) -> pd.DataFrame:
    """Violin distributions of the feasible flux space per reaction."""
    if samples is None or samples.empty:
        return _empty(fig, "No flux-sampling data.")
    # Choose the most variable reactions to display.
    var = samples.var().sort_values(ascending=False)
    cols = list(var.head(top).index)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.violinplot([samples[c].values for c in cols], showmeans=True, vert=False)
    ax.set_yticks(range(1, len(cols) + 1))
    ax.set_yticklabels(cols, fontsize=7)
    ax.axvline(0, color="#B8C0CC", lw=0.8)
    ax.set_xlabel("flux")
    ax.set_title(f"Flux sampling — {len(cols)} most variable reactions")
    theme.style_axes(ax)
    return samples[cols].describe().reset_index(names="stat")


# --- Omics (best-effort; used when expression data is present) ---------------
def volcano(fig, df: pd.DataFrame) -> pd.DataFrame:
    """Volcano plot: log2 fold-change vs −log10 p-value."""
    need = {"log2fc", "pvalue"}
    if df is None or df.empty or not need.issubset(df.columns):
        return _empty(fig, "Volcano needs columns 'log2fc' and 'pvalue'.")
    d = df.copy()
    d["neglogp"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    sig = (d["pvalue"] < 0.05) & (d["log2fc"].abs() > 1)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.scatter(d.loc[~sig, "log2fc"], d.loc[~sig, "neglogp"], s=8, color=theme.NEUTRAL, alpha=0.6)
    ax.scatter(d.loc[sig, "log2fc"], d.loc[sig, "neglogp"], s=10, color=theme.UP)
    ax.axvline(0, color="#B8C0CC", lw=0.8)
    ax.set_xlabel("log₂ fold-change")
    ax.set_ylabel("−log₁₀ p")
    ax.set_title("Volcano plot")
    theme.style_axes(ax)
    return d


def pca_scatter(fig, df: pd.DataFrame) -> pd.DataFrame:
    """PCA of conditions (samples = columns) — first two components."""
    if df is None or df.shape[1] < 2:
        return _empty(fig, "PCA needs a numeric matrix with ≥2 conditions (columns).")
    x = df.select_dtypes("number").dropna()
    if x.shape[1] < 2:
        return _empty(fig, "PCA needs ≥2 numeric conditions.")
    m = x.values.T                      # conditions as rows
    m = m - m.mean(axis=0)
    try:
        u, s, _vt = np.linalg.svd(m, full_matrices=False)
        pcs = u[:, :2] * s[:2]
    except Exception:  # noqa: BLE001
        return _empty(fig, "PCA could not be computed for this data.")
    fig.clear()
    ax = fig.add_subplot(111)
    ax.scatter(pcs[:, 0], pcs[:, 1], s=40, color=theme.CATEGORICAL[0])
    for i, name in enumerate(x.columns):
        ax.annotate(str(name), (pcs[i, 0], pcs[i, 1]), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of conditions")
    theme.style_axes(ax)
    return pd.DataFrame(pcs, index=x.columns, columns=["PC1", "PC2"]).reset_index(names="condition")
