"""Elementary Flux Modes (EFMs) for small subnetworks (Phase 6).

EFMs are the minimal, indivisible steady-state flux routes through a network — the
"building blocks" of metabolism. Enumerating them is exponential, so genome-scale
EFM analysis is intractable; this module therefore operates on a *category /
subnetwork* with a hard reaction-count cap.

Method: split reversible reactions into non-negative forward/backward parts so the
feasible space is the pointed flux cone ``{v >= 0 : S v = 0}``. Its extreme rays,
enumerated by the double-description method (pycddlib), are the elementary flux
modes. Each ray is mapped back to net reaction fluxes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Tuple

import cobra
import pandas as pd

from .. import categories

MAX_REACTIONS = 30  # EFM enumeration is exponential — keep subnetworks small.


class PathwayError(Exception):
    """Raised when EFM enumeration cannot run (too large, missing solver, etc.)."""


@dataclass
class EFMResult:
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    n_modes: int = 0
    note: str = ""


def enumerate_efms(
    model: cobra.Model,
    reaction_ids,
    *,
    max_modes: int = 200,
) -> EFMResult:
    """Enumerate elementary flux modes of the subnetwork defined by ``reaction_ids``.

    A standalone sub-model is built (cut metabolites get free exchanges) and its
    EFMs are computed. Raises :class:`PathwayError` if the subnetwork is too large.
    """
    try:
        import cdd
    except ImportError as exc:  # noqa: BLE001
        raise PathwayError(
            "Elementary flux mode analysis requires the 'pycddlib' package.") from exc

    sub = categories.build_subset_model(model, reaction_ids)
    reactions = list(sub.reactions)
    if len(reactions) > MAX_REACTIONS:
        raise PathwayError(
            f"This subnetwork has {len(reactions)} reactions (incl. boundary exchanges); "
            f"EFM enumeration is limited to {MAX_REACTIONS}. Use a smaller category.")

    metabolites = list(sub.metabolites)
    met_index = {m.id: i for i, m in enumerate(metabolites)}

    # Build split (non-negative) reaction variables.
    split: List[Tuple[str, int]] = []  # (reaction_id, sign)
    for rxn in reactions:
        rev = rxn.lower_bound < 0 < rxn.upper_bound
        if rev:
            split.append((rxn.id, 1))
            split.append((rxn.id, -1))
        elif rxn.upper_bound <= 0 and rxn.lower_bound < 0:
            split.append((rxn.id, -1))
        else:
            split.append((rxn.id, 1))

    n_split = len(split)
    # Stoichiometric matrix columns for split variables.
    # H-representation rows: [b, a...] meaning b + a·v >= 0.
    rows = []
    lin = set()
    # Steady-state equalities: for each metabolite, S·v = 0.
    for m_i, met in enumerate(metabolites):
        row = [Fraction(0)] * (n_split + 1)
        for s_j, (rid, sign) in enumerate(split):
            coeff = sub.reactions.get_by_id(rid).metabolites.get(met, 0)
            if coeff:
                row[1 + s_j] = Fraction(sign) * Fraction(coeff).limit_denominator(10**6)
        rows.append(row)
        lin.add(m_i)
    # Non-negativity inequalities: each split variable >= 0.
    for s_j in range(n_split):
        e = [Fraction(0)] * (n_split + 1)
        e[1 + s_j] = Fraction(1)
        rows.append(e)

    try:
        mat = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY)
        mat.lin_set = lin
        poly = cdd.polyhedron_from_matrix(mat)
        gens = cdd.copy_generators(poly)
    except Exception as exc:  # noqa: BLE001
        raise PathwayError(f"EFM enumeration failed:\n{exc}") from exc

    # Extreme rays (generator rows whose leading entry is 0).
    rays = [row[1:] for row in gens.array if row[0] == 0]
    if not rays:
        return EFMResult(note="No elementary flux modes found for this subnetwork.")

    modes = []
    for ray in rays[:max_modes]:
        net: Dict[str, float] = {}
        for value, (rid, sign) in zip(ray, split):
            if value != 0:
                net[rid] = net.get(rid, 0.0) + sign * float(value)
        net = {k: v for k, v in net.items() if abs(v) > 1e-9}
        if net:
            modes.append(net)

    rows_out = []
    for i, mode in enumerate(modes, start=1):
        rows_out.append({
            "EFM": i,
            "n_reactions": len(mode),
            "reactions": ", ".join(f"{rid}:{v:.3g}" for rid, v in sorted(mode.items())),
        })
    table = pd.DataFrame(rows_out)
    note = f"{len(modes)} elementary flux modes in a {len(reactions)}-reaction subnetwork."
    if len(rays) > max_modes:
        note += f" (showing first {max_modes})"
    return EFMResult(table=table, n_modes=len(modes), note=note)
