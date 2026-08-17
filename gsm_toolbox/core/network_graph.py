"""Build a navigable metabolic network graph from a cobra.Model.

The graph is *bipartite*: metabolite nodes and reaction nodes, with edges
connecting a reaction to each metabolite it consumes/produces. For genome-scale
models a full layout is an unreadable "hairball", so the default mode is a
**focused neighborhood**: pick a seed reaction or metabolite and expand out a
fixed number of hops.

Currency metabolites (ATP, water, protons, NAD(P)H, CO2 ...) appear in hundreds
of reactions and dominate the layout; they are hidden by default.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import cobra
import networkx as nx

# Base names of ubiquitous "currency" metabolites, matched ignoring the trailing
# compartment suffix (e.g. ``atp_c`` -> base ``atp``).
CURRENCY_BASES: Set[str] = {
    "h", "h2o", "atp", "adp", "amp", "pi", "ppi",
    "nad", "nadh", "nadp", "nadph", "co2", "o2",
    "coa", "nh4", "fad", "fadh2", "q8", "q8h2", "accoa",
    # Nucleotide energy/redox carriers
    "gtp", "gdp", "gmp", "ctp", "cdp", "cmp", "utp", "udp", "ump",
    "itp", "idp", "imp", "dttp", "datp", "dctp", "dgtp",
    # Redox / one-carbon / group-transfer cofactors that decorate reactions but
    # are not the carbon skeleton being transformed (so they belong beside the
    # arrow, not as a principal metabolite in a pathway drawing).
    "fdxrd", "fdxox", "fdxo", "fdxr", "fdxh", "fldrd", "fldox",   # ferredoxin/flavodoxin
    "trdrd", "trdox",                                            # thioredoxin
    "gthrd", "gthox",                                            # glutathione
    "cbl1", "cbl2", "adocbl", "aqcbl", "b12",                    # cobalamin (B12)
    "amet", "ahcys",                                            # SAM / SAH
    "acp",                                                       # acyl-carrier protein
    "thf", "5mthf", "methf", "mlthf", "10fthf", "5fthf",         # folates
    "mqn8", "mql8", "2dmmql8", "q8h2",                           # (mena)quinones
    "so4", "so3", "h2s", "hco3", "h2o2", "fmn", "fmnh2",
    "nadphx", "biotin", "btn", "tpp", "thmpp", "pydx5p",         # biotin, thiamine-PP, PLP
}

# Short, familiar display names for common metabolites, shown instead of long
# database names (e.g. "Nicotinamide adenine dinucleotide - reduced" -> "NADH").
SHORT_NAMES = {
    "h": "H+", "h2o": "H2O", "atp": "ATP", "adp": "ADP", "amp": "AMP",
    "pi": "Pi", "ppi": "PPi", "nad": "NAD+", "nadh": "NADH", "nadp": "NADP+",
    "nadph": "NADPH", "co2": "CO2", "o2": "O2", "coa": "CoA", "nh4": "NH4+",
    "fad": "FAD", "fadh2": "FADH2", "q8": "Q8", "q8h2": "Q8H2",
    "accoa": "Acetyl-CoA", "gtp": "GTP", "gdp": "GDP", "utp": "UTP", "udp": "UDP",
    "ctp": "CTP", "cdp": "CDP", "fmn": "FMN", "fmnh2": "FMNH2", "so4": "Sulfate",
    "h2o2": "H2O2", "hco3": "HCO3-", "nadphx": "NADPHX",
    "fdxrd": "reduced ferredoxin", "fdxox": "oxidized ferredoxin",
    "fdxo": "oxidized ferredoxin", "fdxr": "reduced ferredoxin",
    "trdrd": "reduced thioredoxin", "trdox": "oxidized thioredoxin",
    "gthrd": "glutathione (red)", "gthox": "glutathione (ox)",
    "amet": "SAM", "ahcys": "SAH", "acp": "ACP", "thf": "THF",
    "adocbl": "adenosylcobalamin", "cbl1": "cob(I)alamin", "cbl2": "cob(II)alamin",
    "so3": "sulfite", "h2s": "H2S", "biotin": "biotin", "btn": "biotin",
}


def clean_label(text: str) -> str:
    """Tidy a reaction/metabolite label for display: decode the SBML-safe
    ``_LPAREN_e_RPAREN_`` compartment encoding (old BiGG style) back to ``(e)`` and
    drop the leftover ``LPAREN``/``RPAREN`` tokens that otherwise leak into the UI
    (e.g. iJN678's ``EX hco3 LPAREN e RPAREN``). Pure display cleanup — never used
    to change a reaction/metabolite ``id`` or identity."""
    s = str(text or "")
    s = s.replace("_LPAREN_", "(").replace("_RPAREN_", ")")
    s = re.sub(r"[ _]*LPAREN[ _]*", "(", s)
    s = re.sub(r"[ _]*RPAREN[ _]*", ")", s)
    s = s.replace("( ", "(").replace(" )", ")")
    return re.sub(r"\s+", " ", s).strip()


def short_metabolite_name(met_id: str, name: str = "") -> str:
    """A concise display name: a short standard name for common metabolites,
    otherwise the given name (falling back to the id). The result is tidied of the
    SBML LPAREN/RPAREN compartment encoding so it is always UI-safe."""
    base = met_id.rsplit("_", 1)[0].lower() if "_" in met_id else met_id.lower()
    if base in SHORT_NAMES:
        return SHORT_NAMES[base]
    return clean_label(name or met_id)


def display_reaction_name(rxn) -> str:
    """The display name of a reaction, tidied of SBML LPAREN/RPAREN encoding."""
    return clean_label(getattr(rxn, "name", "") or getattr(rxn, "id", "") or "")


def reaction_equation(rxn) -> str:
    """A human-readable reaction equation using short/familiar metabolite names
    (ATP, NADH…) for currency metabolites and names for the rest."""
    reactants, products = [], []
    for met, coeff in rxn.metabolites.items():
        label = short_metabolite_name(met.id, met.name or "")
        term = label if abs(coeff) == 1 else f"{abs(coeff):g} {label}"
        (reactants if coeff < 0 else products).append(term)
    arrow = "<=>" if rxn.reversibility else "-->"
    return f"{' + '.join(reactants)} {arrow} {' + '.join(products)}"


@dataclass
class GraphNode:
    node_id: str
    label: str
    kind: str  # "metabolite" | "reaction"
    x: float = 0.0
    y: float = 0.0
    data: dict = field(default_factory=dict)


@dataclass
class GraphEdge:
    source: str  # node_id
    target: str  # node_id
    stoichiometry: float = 1.0
    flux: Optional[float] = None


@dataclass
class NetworkGraph:
    nodes: List[GraphNode]
    edges: List[GraphEdge]

    def node_count(self) -> int:
        return len(self.nodes)


def _met_base(met_id: str) -> str:
    """Strip a trailing ``_<compartment>`` suffix to get the metabolite base name."""
    return met_id.rsplit("_", 1)[0] if "_" in met_id else met_id


def _is_currency(met: cobra.Metabolite) -> bool:
    return _met_base(met.id).lower() in CURRENCY_BASES


def reaction_type(rxn: cobra.Reaction) -> str:
    """Classify a reaction for iconography: exchange | transport | reversible | irreversible."""
    if rxn.boundary:
        return "exchange"
    # Transport: the SAME metabolite (same base id) appears in >1 compartment —
    # i.e. a species is moved across a membrane. Merely involving metabolites from
    # different compartments (e.g. a hydrolase written with a cross-compartment
    # cofactor) is NOT transport, so we don't use the coarse "compartments > 1" test.
    comps_by_base: Dict[str, Set[str]] = {}
    for met in rxn.metabolites:
        comps_by_base.setdefault(_met_base(met.id), set()).add(met.compartment or "")
    if any(len(comps) > 1 for comps in comps_by_base.values()):
        return "transport"
    return "reversible" if rxn.reversibility else "irreversible"


def build_graph(
    model: cobra.Model,
    *,
    center: Optional[str] = None,
    radius: int = 1,
    hide_currency: bool = True,
    max_nodes: int = 400,
    fluxes: Optional[Dict[str, float]] = None,
    reaction_subset: Optional[Set[str]] = None,
    seed_reactions: Optional[Set[str]] = None,
    layout: str = "radial",
    layout_seed: int = 42,
) -> NetworkGraph:
    """Build a (optionally neighborhood-limited) network graph with a 2-D layout.

    Parameters
    ----------
    center:
        A reaction or metabolite ID to center the view on. If ``None`` the whole
        model is used (capped at ``max_nodes`` reactions to stay responsive).
    radius:
        Number of hops to expand around ``center`` (only used when ``center`` set).
    hide_currency:
        Drop ubiquitous currency metabolites to keep the map readable.
    fluxes:
        Optional ``{reaction_id: flux}`` mapping (e.g. from an FBA result) used to
        annotate edges for overlay rendering.
    reaction_subset:
        If given, the graph is built from exactly these reactions (used to view a
        category/pathway in isolation), ignoring ``center``/``radius``.
    """
    g = nx.Graph()
    if seed_reactions is not None:
        # Expand outward from a set of seed reactions (e.g. a category) over the
        # whole network by ``radius`` hops; radius 0 keeps only the seeds.
        reactions = _expand_from_reactions(
            model, seed_reactions, radius, hide_currency, max_nodes)
    elif reaction_subset is not None:
        reactions = [model.reactions.get_by_id(r) for r in reaction_subset
                     if model.reactions.has_id(r)]
    else:
        reactions = _select_reactions(model, center, radius, hide_currency, max_nodes)

    for rxn in reactions:
        r_node = f"R:{rxn.id}"
        g.add_node(r_node, label=rxn.id, kind="reaction", name=rxn.name or "",
                   rxn_type=reaction_type(rxn))
        for met, coeff in rxn.metabolites.items():
            if hide_currency and _is_currency(met):
                continue
            m_node = f"M:{met.id}"
            if m_node not in g:
                # ``label`` stays the id (used for identity/click-resolution);
                # ``display`` is what the map shows — currency metabolites get their
                # short familiar name (ATP, NADH…), others show their id.
                g.add_node(m_node, label=met.id, kind="metabolite", name=met.name or "",
                           display=short_metabolite_name(met.id, ""),
                           compartment=met.compartment or "")
            flux = None if fluxes is None else fluxes.get(rxn.id)
            g.add_edge(r_node, m_node, stoichiometry=float(coeff), flux=flux)

    # Drop metabolite nodes that ended up isolated (e.g. all neighbors hidden).
    isolated = [n for n in g.nodes if g.degree(n) == 0]
    g.remove_nodes_from(isolated)

    center_key = None
    if center is not None:
        if model.reactions.has_id(center):
            center_key = f"R:{center}"
        elif model.metabolites.has_id(center):
            center_key = f"M:{center}"
    pos = _layout(g, mode=layout, center_key=center_key, seed=layout_seed)

    nodes = [
        GraphNode(
            node_id=n,
            label=attrs.get("label", n),
            kind=attrs.get("kind", "metabolite"),
            x=float(pos[n][0]),
            y=float(pos[n][1]),
            data={"name": attrs.get("name", ""), "rxn_type": attrs.get("rxn_type", ""),
                  "display": attrs.get("display", attrs.get("label", n)),
                  "compartment": attrs.get("compartment", "")},
        )
        for n, attrs in g.nodes(data=True)
    ]
    edges = [
        GraphEdge(
            source=u,
            target=v,
            stoichiometry=attrs.get("stoichiometry", 1.0),
            flux=attrs.get("flux"),
        )
        for u, v, attrs in g.edges(data=True)
    ]
    return NetworkGraph(nodes=nodes, edges=edges)


def _select_reactions(
    model: cobra.Model,
    center: Optional[str],
    radius: int,
    hide_currency: bool,
    max_nodes: int,
) -> List[cobra.Reaction]:
    if center is None:
        # Whole-model view, capped for responsiveness.
        return list(model.reactions)[:max_nodes]

    # Resolve the seed (reaction or metabolite) and BFS outward by `radius` hops.
    seed_reactions: Set[str] = set()
    frontier_mets: Set[str] = set()

    if model.reactions.has_id(center):
        seed_reactions.add(center)
        for met in model.reactions.get_by_id(center).metabolites:
            if not (hide_currency and _is_currency(met)):
                frontier_mets.add(met.id)
    elif model.metabolites.has_id(center):
        frontier_mets.add(center)
    else:
        raise KeyError(f"'{center}' is neither a reaction nor a metabolite ID.")

    visited_mets: Set[str] = set()
    for _ in range(max(radius, 1)):
        next_mets: Set[str] = set()
        for met_id in frontier_mets:
            if met_id in visited_mets:
                continue
            visited_mets.add(met_id)
            met = model.metabolites.get_by_id(met_id)
            for rxn in met.reactions:
                seed_reactions.add(rxn.id)
                for m in rxn.metabolites:
                    if hide_currency and _is_currency(m):
                        continue
                    if m.id not in visited_mets:
                        next_mets.add(m.id)
        frontier_mets = next_mets

    return [model.reactions.get_by_id(r) for r in seed_reactions][:max_nodes]


def _expand_from_reactions(
    model: cobra.Model,
    seeds: Set[str],
    radius: int,
    hide_currency: bool,
    max_nodes: int,
) -> List[cobra.Reaction]:
    """Start from ``seeds`` reactions and add reactions within ``radius`` hops over
    the whole model (radius 0 keeps just the seeds)."""
    chosen: Set[str] = {r for r in seeds if model.reactions.has_id(r)}
    frontier_mets: Set[str] = set()
    for rid in chosen:
        for met in model.reactions.get_by_id(rid).metabolites:
            if not (hide_currency and _is_currency(met)):
                frontier_mets.add(met.id)

    visited_mets: Set[str] = set()
    for _ in range(max(radius, 0)):
        next_mets: Set[str] = set()
        for met_id in frontier_mets:
            if met_id in visited_mets:
                continue
            visited_mets.add(met_id)
            met = model.metabolites.get_by_id(met_id)
            for rxn in met.reactions:
                chosen.add(rxn.id)
                for m in rxn.metabolites:
                    if hide_currency and _is_currency(m):
                        continue
                    if m.id not in visited_mets:
                        next_mets.add(m.id)
            if len(chosen) >= max_nodes:
                break
        frontier_mets = next_mets
        if len(chosen) >= max_nodes:
            break
    return [model.reactions.get_by_id(r) for r in chosen][:max_nodes]


def _layout(g: "nx.Graph", *, mode: str, center_key: Optional[str], seed: int) -> Dict[str, tuple]:
    """Compute 2-D node positions.

    * ``radial``  — concentric shells by hop-distance from the focus node; intuitive
      for exploring the neighborhood of one reaction/metabolite.
    * ``layered`` — left-to-right layers by hop-distance (a flow-like view that
      reads like a pathway diagram).
    * ``force``   — force-directed spring layout (good for whole-model overviews).
    """
    n = g.number_of_nodes()
    if n == 0:
        return {}

    if mode == "radial" and center_key in g:
        try:
            lengths = nx.single_source_shortest_path_length(g, center_key)
            max_d = max(lengths.values()) if lengths else 0
            shells = [[] for _ in range(max_d + 1)]
            for node in g.nodes:
                shells[lengths.get(node, max_d)].append(node)
            shells = [s for s in shells if s]
            return nx.shell_layout(g, nlist=shells)
        except Exception:  # noqa: BLE001
            pass

    if mode == "layered":
        try:
            root = center_key if center_key in g else max(g.degree, key=lambda kv: kv[1])[0]
            lengths = nx.single_source_shortest_path_length(g, root)
            for node in g.nodes:
                g.nodes[node]["layer"] = lengths.get(node, 0)
            return nx.multipartite_layout(g, subset_key="layer", align="vertical")
        except Exception:  # noqa: BLE001
            pass

    # Force-directed: spread nodes out more for readability on larger graphs.
    try:
        return nx.spring_layout(g, seed=seed, k=2.5 / (n ** 0.5), iterations=80)
    except Exception:  # noqa: BLE001
        return nx.circular_layout(g)
