"""One theming module shared by every figure (the proposal's cross-cutting item).

Colour-blind-safe defaults — ``viridis`` for magnitude (sequential) and ``RdBu_r``
for differences (diverging) — plus a consistent categorical palette, typography
and axis styling so all plots read as one system in light or dark contexts.
"""

from __future__ import annotations

# Colour-blind-safe colormaps (matplotlib names).
SEQUENTIAL = "viridis"          # magnitude, e.g. |flux|
DIVERGING = "RdBu_r"            # differences, e.g. Δflux (red = more, blue = less)

# A colour-blind-friendly categorical palette (Okabe–Ito).
CATEGORICAL = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]

# Semantic colours matching the network map's flux convention.
UP = "#D55E00"       # more flux / amplification (warm)
DOWN = "#0072B2"     # less flux / knock-down (cool)
NEUTRAL = "#9AA6B2"


def apply_style(fig=None) -> None:
    """Apply consistent typography/grid styling to a matplotlib Figure (or the
    current rcParams if ``fig`` is None). Safe to call repeatedly."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "600",
        "axes.labelsize": 9,
        "axes.edgecolor": "#B8C0CC",
        "axes.grid": True,
        "grid.color": "#E3E8EF",
        "grid.linewidth": 0.6,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "figure.autolayout": True,
        "svg.fonttype": "none",     # keep text editable in exported SVG/PDF
    })
    if fig is not None:
        fig.set_facecolor("white")


def style_axes(ax) -> None:
    """De-clutter an Axes: hide top/right spines, soften the grid."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, alpha=0.35)
    ax.tick_params(labelsize=8)
