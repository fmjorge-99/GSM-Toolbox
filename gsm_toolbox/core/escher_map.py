"""Generate Escher-compatible metabolic **map JSON** from any cobra model.

Escher (https://escher.github.io) is the field-standard interactive metabolic-map
renderer, but it needs a *pre-drawn map* — a JSON describing where every reaction
arrow and metabolite circle sits. Curated maps only exist for a handful of BiGG
models, so to let the ToolBox draw *any* loaded model (including the user's own
non-standard / cyanobacterial models) we synthesise a map automatically: we reuse
the same graph + 2-D layout the native Network Map uses (:mod:`network_graph`) and
translate it into Escher's node/segment schema.

The output is the two-element ``[metadata, body]`` list Escher's ``Builder``
expects, with ``reactions``/``nodes``/``segments``/``canvas`` populated. Flux is
*not* baked into the map — the panel feeds it to the Builder as ``reaction_data``
at runtime, which is what makes the overlay/difference/FBA-navigator modes live.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Set

import cobra
import networkx as nx

from . import network_graph

# Escher's schema URL for the 1-0-0 map format (what the bundled renderer reads).
_SCHEMA = "https://escher.github.io/escher/jsonschema/1-0-0#"

# Layout scaling: the graph layout returns coordinates in roughly [-1, 1]; Escher
# works in a large pixel space where reaction arrows are a few hundred px long. A
# generous spacing keeps connectors long and leaves room for the (enlarged) labels.
_NODE_SPACING = 340.0
_MET_LABEL_DX = 16.0
_MET_LABEL_DY = -22.0
_RXN_LABEL_DY = -36.0


def _display_name(node: network_graph.GraphNode) -> str:
    disp = node.data.get("display") or node.data.get("name") or node.label
    return str(disp)


def _clean_layout(edges, node_ids, seed: int = 42) -> Dict[str, tuple]:
    """Compute node positions that spread the network out with as few crossings as
    practical. Kamada–Kawai gives the tidiest, most 'geometric' arrangement for the
    small/medium subnetworks we recommend focusing on; for large whole-model graphs
    we fall back to a well-separated force layout (higher k, more iterations)."""
    g = nx.Graph()
    g.add_nodes_from(node_ids)
    for u, v in edges:
        g.add_edge(u, v)
    n = g.number_of_nodes()
    if n == 0:
        return {}
    if n <= 200:
        try:
            # Kamada–Kawai is deterministic and untangles medium graphs nicely; seed
            # it with a circular layout so disconnected pieces don't pile up.
            return nx.kamada_kawai_layout(g, pos=nx.circular_layout(g))
        except Exception:  # noqa: BLE001 - fall back if scipy path fails
            pass
    try:
        return nx.spring_layout(g, seed=seed, k=3.2 / (n ** 0.5), iterations=200)
    except Exception:  # noqa: BLE001
        return nx.circular_layout(g)


def build_escher_map(
    model: cobra.Model,
    reaction_ids: Optional[Iterable[str]] = None,
    *,
    name: str = "Auto-generated map",
    hide_currency: bool = True,
    center: Optional[str] = None,
    radius: int = 1,
    layout: str = "force",
    max_nodes: int = 400,
    layout_seed: int = 42,
) -> list:
    """Build an Escher map JSON (``[metadata, body]``) for the given reactions.

    Parameters mirror :func:`network_graph.build_graph`. If ``reaction_ids`` is
    given the map is restricted to exactly those reactions (a category/pathway);
    otherwise ``center``/``radius`` (or the whole model, capped) are used.
    """
    subset: Optional[Set[str]] = set(reaction_ids) if reaction_ids else None
    ng = network_graph.build_graph(
        model,
        center=center,
        radius=radius,
        hide_currency=hide_currency,
        max_nodes=max_nodes,
        reaction_subset=subset,
        layout=layout,
        layout_seed=layout_seed,
    )

    # --- compute a clean layout, then scale into Escher pixel space ----------
    raw = _clean_layout([(e.source, e.target) for e in ng.edges],
                        [n.node_id for n in ng.nodes], seed=layout_seed)
    if ng.nodes and raw:
        xs = [raw[n.node_id][0] for n in ng.nodes if n.node_id in raw]
        ys = [raw[n.node_id][1] for n in ng.nodes if n.node_id in raw]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        span = max(max_x - min_x, max_y - min_y) or 1.0
        # spread nodes ~_NODE_SPACING apart on average
        scale = _NODE_SPACING * math.sqrt(max(len(ng.nodes), 1)) / span
    else:
        min_x = min_y = 0.0
        scale = 1.0

    def _px(node_id: str):
        x, y = raw.get(node_id, (0.0, 0.0))
        return (round((x - min_x) * scale + 100.0, 2),
                round((y - min_y) * scale + 100.0, 2))

    pos: Dict[str, tuple] = {n.node_id: _px(n.node_id) for n in ng.nodes}

    # --- allocate stable string uids -----------------------------------------
    uid = [0]

    def _next() -> str:
        uid[0] += 1
        return str(uid[0])

    nodes: Dict[str, dict] = {}
    reactions: Dict[str, dict] = {}

    met_node_uid: Dict[str, str] = {}   # graph "M:<id>" -> escher node uid
    met_id_of: Dict[str, str] = {}      # graph node_id -> metabolite id

    # metabolite nodes
    for n in ng.nodes:
        if n.kind != "metabolite":
            continue
        x, y = pos[n.node_id]
        u = _next()
        met_node_uid[n.node_id] = u
        met_id_of[n.node_id] = n.label
        nodes[u] = {
            "node_type": "metabolite",
            "x": x, "y": y,
            "bigg_id": n.label,
            "name": _display_name(n),
            "label_x": x + _MET_LABEL_DX,
            "label_y": y + _MET_LABEL_DY,
            "node_is_primary": True,
        }

    # adjacency: reaction graph-node -> list of connected metabolite graph-nodes
    adj: Dict[str, List[str]] = {}
    for e in ng.edges:
        r, m = (e.source, e.target) if e.source.startswith("R:") else (e.target, e.source)
        adj.setdefault(r, []).append(m)

    for n in ng.nodes:
        if n.kind != "reaction":
            continue
        rid = n.label
        if not model.reactions.has_id(rid):
            continue
        rxn = model.reactions.get_by_id(rid)
        mx, my = pos[n.node_id]
        mid_uid = _next()
        nodes[mid_uid] = {"node_type": "midmarker", "x": mx, "y": my}

        coeff_by_met = {met.id: float(c) for met, c in rxn.metabolites.items()}
        segments: Dict[str, dict] = {}
        for m_graph in adj.get(n.node_id, []):
            m_uid = met_node_uid.get(m_graph)
            if m_uid is None:
                continue
            met_id = met_id_of[m_graph]
            coeff = coeff_by_met.get(met_id, -1.0)
            mxn, myn = pos[m_graph]
            # reactant (coeff<0): metabolite -> midmarker; product: midmarker -> met,
            # so the arrowhead lands on the product side.
            if coeff < 0:
                frm, to = m_uid, mid_uid
                ax, ay, bx, by = mxn, myn, mx, my
            else:
                frm, to = mid_uid, m_uid
                ax, ay, bx, by = mx, my, mxn, myn
            b1 = {"x": round(ax + (bx - ax) / 3.0, 2), "y": round(ay + (by - ay) / 3.0, 2)}
            b2 = {"x": round(ax + 2 * (bx - ax) / 3.0, 2), "y": round(ay + 2 * (by - ay) / 3.0, 2)}
            segments[_next()] = {"from_node_id": frm, "to_node_id": to, "b1": b1, "b2": b2}

        # Place the reaction label just above its midmarker (a white halo in the
        # renderer keeps it legible where it crosses a connector).
        genes = [{"bigg_id": g.id, "name": g.name or g.id} for g in rxn.genes]
        reactions[_next()] = {
            "name": rxn.name or rid,
            "bigg_id": rid,
            "reversibility": bool(rxn.reversibility),
            "label_x": mx + 6.0,
            "label_y": my + _RXN_LABEL_DY,
            "gene_reaction_rule": rxn.gene_reaction_rule or "",
            "genes": genes,
            "metabolites": [
                {"coefficient": float(c), "bigg_id": met.id}
                for met, c in rxn.metabolites.items()
            ],
            "segments": segments,
        }

    # --- canvas bounding box (with margin) -----------------------------------
    if nodes:
        allx = [nd["x"] for nd in nodes.values()]
        ally = [nd["y"] for nd in nodes.values()]
        margin = 200.0
        canvas = {
            "x": min(allx) - margin,
            "y": min(ally) - margin,
            "width": (max(allx) - min(allx)) + 2 * margin,
            "height": (max(ally) - min(ally)) + 2 * margin,
        }
    else:
        canvas = {"x": 0, "y": 0, "width": 1000, "height": 1000}

    meta = {
        "map_name": name,
        "map_id": f"gsmtb-{abs(hash(name)) & 0xffffffff:08x}",
        "map_description": "Auto-generated by GSM ToolBox",
        "homepage": "https://escher.github.io",
        "schema": _SCHEMA,
    }
    body = {"reactions": reactions, "nodes": nodes,
            "text_labels": {}, "canvas": canvas}
    return [meta, body]


def flux_overlay(fluxes: Dict[str, float]) -> Dict[str, float]:
    """Coerce a flux mapping into the ``{bigg_id: value}`` dict Escher's
    ``reaction_data`` expects (drops NaN/None, keeps sign for direction)."""
    out: Dict[str, float] = {}
    for rid, v in (fluxes or {}).items():
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if fv != fv:  # NaN
            continue
        out[str(rid)] = fv
    return out
