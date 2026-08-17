"""Flux sampling — the feasible flux space as distributions.

Samples the solution space (COBRApy's OptGP sampler) to convey uncertainty and
alternate optima: violin/ridgeline distributions per reaction (see the graphical
engine's per-utility catalogue).
"""

from __future__ import annotations

from typing import List, Optional

import cobra
import pandas as pd


class SamplingError(Exception):
    """Raised when flux sampling cannot run."""


def run_flux_sampling(model: cobra.Model, *, n_samples: int = 500,
                      reaction_ids: Optional[List[str]] = None,
                      thinning: int = 100) -> pd.DataFrame:
    """Return an ``n_samples × reactions`` DataFrame of sampled fluxes.

    ``reaction_ids`` restricts the returned columns (e.g. a category) so the
    result stays legible; sampling itself always runs on the whole model. Uses
    OptGP (parallel-friendly) with a single process for determinism on Windows.
    """
    from cobra.sampling import sample

    try:
        samples = sample(model, n=int(n_samples), method="optgp",
                         thinning=int(thinning), processes=1)
    except Exception as exc:  # noqa: BLE001
        raise SamplingError(
            f"Flux sampling failed: {exc}. The model may be infeasible under the "
            "current medium, or the range is degenerate.") from exc

    if reaction_ids:
        keep = [r for r in reaction_ids if r in samples.columns]
        if keep:
            samples = samples[keep]
    return samples
