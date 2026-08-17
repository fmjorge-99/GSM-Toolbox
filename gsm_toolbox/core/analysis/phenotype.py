"""Phenotype / robustness analyses: production envelope, robustness scan, PhPP.

These explore how the objective (usually growth) and product yields respond to
changing one or two reaction fluxes — the workhorse plots of metabolic
engineering yield analysis.
"""

from __future__ import annotations

from typing import List, Optional

import cobra
import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, production_envelope


def run_production_envelope(
    model: cobra.Model,
    target_reaction: str,
    *,
    objective: Optional[str] = None,
    points: int = 20,
) -> pd.DataFrame:
    """Production envelope: feasible range of ``target_reaction`` flux vs growth.

    Returns the DataFrame from cobra's ``production_envelope`` (columns include the
    target flux, and the min/max of the objective at each level).
    """
    return production_envelope(
        model, [target_reaction], objective=objective, points=points
    )


def run_robustness(
    model: cobra.Model,
    control_reaction: str,
    *,
    points: int = 25,
    lower: Optional[float] = None,
    upper: Optional[float] = None,
) -> pd.DataFrame:
    """Robustness analysis: fix ``control_reaction`` across a range, record objective.

    Reveals "tipping points" where forcing a reaction's flux crashes growth.
    """
    rxn = model.reactions.get_by_id(control_reaction)
    # Default scan range: the reaction's feasible flux span (FVA), so the scan
    # stays inside feasible territory instead of hitting infeasible extremes
    # (which previously produced a near-empty plot).
    if lower is None or upper is None:
        fva = flux_variability_analysis(model, [rxn], fraction_of_optimum=0.0)
        fva_lo, fva_hi = float(fva["minimum"].iloc[0]), float(fva["maximum"].iloc[0])
        lo = fva_lo if lower is None else lower
        hi = fva_hi if upper is None else upper
    else:
        lo, hi = lower, upper
    if not np.isfinite(lo):
        lo = -1000.0
    if not np.isfinite(hi):
        hi = 1000.0
    if abs(hi - lo) < 1e-9:
        lo, hi = lo - 1.0, hi + 1.0
    values = np.linspace(lo, hi, points)

    rows = []
    for v in values:
        with model:
            rxn.bounds = (v, v)
            sol = model.slim_optimize()
            rows.append({"control_flux": float(v),
                         "objective": float(sol) if (sol is not None and np.isfinite(sol)) else np.nan})
    return pd.DataFrame(rows)


def run_phase_plane(
    model: cobra.Model,
    reaction_x: str,
    reaction_y: str,
    *,
    points: int = 15,
) -> pd.DataFrame:
    """Phenotypic Phase Plane: vary two reactions, record the objective surface.

    Returns a long-format DataFrame with columns ``x``, ``y``, ``objective`` so the
    GUI can render a heatmap / 3-D surface.
    """
    rx = model.reactions.get_by_id(reaction_x)
    ry = model.reactions.get_by_id(reaction_y)
    xs = np.linspace(_finite(rx.lower_bound, -1000), _finite(rx.upper_bound, 1000), points)
    ys = np.linspace(_finite(ry.lower_bound, -1000), _finite(ry.upper_bound, 1000), points)

    rows = []
    for xv in xs:
        for yv in ys:
            with model:
                rx.bounds = (xv, xv)
                ry.bounds = (yv, yv)
                obj = model.slim_optimize()
                rows.append({"x": xv, "y": yv, "objective": float(obj) if obj is not None else np.nan})
    return pd.DataFrame(rows)


def _finite(value: float, fallback: float) -> float:
    return value if np.isfinite(value) else fallback
