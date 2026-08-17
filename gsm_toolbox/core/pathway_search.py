"""Rebuilt heterologous-pathway search engine.

The previous engine only matched a database to the host by exact identifiers, so
native precursors (DMAP, FPP, …) that carried a different id in the database were
not recognised as available — forcing long de-novo routes (e.g. the whole
mevalonate pathway for isoprene instead of the single native DMAP → isoprene
step). It also could not find longer targets end-to-end.

This module fixes that with three pieces:

1. **Compound identity** (:class:`Canonicalizer`) — a union-find over *every*
   cross-reference (InChIKey, KEGG, BiGG, MetaNetX, SEED, ChEBI, and the
   compartment-stripped base id) of every metabolite in the host and the
   database(s). Two metabolites that share *any* reference are the *same
   compound*, regardless of namespace or compartment. So the host's native
   metabolites and the database's are unified.

2. **Two search algorithms** over that identity space (the user chooses):
   * ``retro`` — retrosynthetic AND/OR search from the target back to native
     compounds; returns the *minimal* set of heterologous reactions that connects
     the target to what the host already makes (so a native precursor ends a
     branch immediately).
   * ``expansion`` — forward network expansion (scope) from the host's compounds.

3. **Database merging** (:func:`merge_databases`) — dedup metabolites (by identity,
   collapsed to a single cytosolic copy) and reactions (by stoichiometry) into one
   unified in-memory database, so a compound/reaction exists exactly once.

Currency metabolites (ATP, H2O, NAD(P)H, CoA, …) are always treated as available.
Everything is compartment-agnostic during the search.
"""

from __future__ import annotations

import os
import weakref
from typing import Dict, List, Optional, Set, Tuple

import cobra

from . import namespace
from .network_graph import CURRENCY_BASES


def _base(mid: str) -> str:
    return mid.rsplit("_", 1)[0] if "_" in mid else mid


def _is_currency_base(mid: str) -> bool:
    return _base(mid).lower() in CURRENCY_BASES


# --- compound identity (union-find over cross-references) --------------------
class Canonicalizer:
    """Assigns a canonical key to a metabolite, unifying compounds that share any
    cross-reference across the host and database(s), ignoring compartments."""

    def __init__(self, models: List[cobra.Model]):
        self._parent: Dict[str, str] = {}
        # key() is called once per metabolite PER REACTION — on a merged database that
        # is ~300k calls, each re-deriving the cross-reference tokens from scratch.
        # Memoise per metabolite object. `models` is retained so those objects stay
        # alive, which keeps id() stable and the cache sound.
        self._models = list(models)
        self._key_cache: Dict[int, str] = {}
        self._tok_cache: Dict[int, Set[str]] = {}
        for model in models:
            for met in model.metabolites:
                toks = self._tokens(met)
                it = iter(toks)
                try:
                    first = next(it)
                except StopIteration:
                    continue
                for t in it:
                    self._union(first, t)

    def _tokens(self, met: cobra.Metabolite) -> Set[str]:
        cached = self._tok_cache.get(id(met))
        if cached is not None:
            return cached
        toks = set(namespace.metabolite_tokens(met))
        toks.add("base:" + _base(met.id).lower())    # same base id also unifies
        self._tok_cache[id(met)] = toks
        return toks

    def _find(self, x: str) -> str:
        self._parent.setdefault(x, x)
        root = x
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[x] != root:      # path compression
            self._parent[x], x = root, self._parent[x]
        return root

    def _union(self, a: str, b: str) -> None:
        ra, rb = self._find(a), self._find(b)
        if ra != rb:
            # keep the lexicographically smaller root for determinism
            lo, hi = sorted((ra, rb))
            self._parent[hi] = lo

    def key(self, met: cobra.Metabolite) -> str:
        cached = self._key_cache.get(id(met))
        if cached is not None:
            return cached
        roots = {self._find(t) for t in self._tokens(met)}
        k = min(roots) if roots else "base:" + _base(met.id).lower()
        self._key_cache[id(met)] = k
        return k


def _currency_keys(canon: Canonicalizer, models: List[cobra.Model]) -> Set[str]:
    """Canonical keys of every currency/cofactor compound across the given models.

    Detection is namespace-agnostic (id base OR metabolite name): a BiGG-style id test
    alone only recognises `fad`/`nadh`/`h2o`-style ids, so ModelSEED's "Reduced flavin"
    or MetaNetX's "MNXM13" looked like ordinary carbon skeletons that had to be
    synthesised from scratch. Any reaction consuming one then never fired, which is why
    e.g. valencene -> nootkatol -> nootkatone (flavin-monooxygenase steps) was
    unreachable even with the reactions present in the database.
    """
    from .pathway_design import _is_currency_met

    keys = set()
    for model in models:
        for met in model.metabolites:
            if _is_currency_met(met):
                keys.add(canon.key(met))
    return keys


def _resolve_target(host: cobra.Model, db: cobra.Model, target_id: str):
    """Find the target metabolite object in the host or database (any compartment)."""
    for model in (host, db):
        if model.metabolites.has_id(target_id):
            return model.metabolites.get_by_id(target_id)
    tb = _base(target_id).lower()
    for model in (host, db):
        for met in model.metabolites:
            if _base(met.id).lower() == tb:
                return met
    return None


# --- retrosynthetic search ---------------------------------------------------
# Cache of the host's producible-compound set so the FVA below is paid once per model
# rather than once per search (and per alternative route). Keyed on the model OBJECT
# via a weak map — keying on id() would be unsound, since a freed model's address can
# be reused by a different one. The stored shape guards against a stale hit after the
# model is edited.
_PRODUCIBLE_CACHE: "weakref.WeakKeyDictionary" = weakref.WeakKeyDictionary()


def _model_fingerprint(host: cobra.Model) -> str:
    """A stable hash of the FVA-relevant model state (reactions, bounds, objective).

    Two models with the same reactions/bounds/objective have the same producible set,
    so a disk cache keyed on this fingerprint is paid once per model state, ever."""
    import hashlib
    h = hashlib.md5()
    for r in sorted(host.reactions, key=lambda x: x.id):
        h.update(f"{r.id}|{r.lower_bound}|{r.upper_bound}|".encode())
    try:
        h.update(str(host.objective.expression).encode())
    except Exception:  # noqa: BLE001
        pass
    return h.hexdigest()


def _producible_cache_path(fp: str) -> str:
    from .cache import databases_dir
    d = os.path.join(os.path.dirname(databases_dir()), "producible")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{fp}.json")


def producible_keys(host: cobra.Model, canon: "Canonicalizer",
                    report: Optional[dict] = None) -> Optional[Set[str]]:
    """Canonical keys of the compounds the host can ACTUALLY make, or ``None`` if this
    could not be determined (caller then treats every present compound as available).

    Being present in the model is not the same as being makeable: iJN678 contains
    4-hydroxybenzyl alcohol whose only producers carry no flux, so its max production is
    zero. This is a *diagnostic* / second-pass check — it must NEVER be used to hide a
    candidate route (a route grounding on a non-producible precursor is still shown, just
    flagged as likely-zero-flux). ``report`` receives ``status`` ∈ {"ok", "fva_failed",
    "no_evidence", "no_reactions"} so the caller can WARN visibly instead of silently
    degrading. Results are cached on disk by a model fingerprint, so the FVA is paid at
    most once per model state.
    """
    if report is None:
        report = {}
    if not host.reactions:
        report["status"] = "no_reactions"
        return None
    fp = _model_fingerprint(host)
    cached = _PRODUCIBLE_CACHE.get(host)
    if cached is not None and cached[0] == fp:
        report["status"] = "ok"
        return cached[1]
    # Disk cache: producible base-keys are namespace-independent, so store them and
    # re-derive canonical keys against THIS canonicalizer on load.
    disk = _producible_cache_path(fp)
    if os.path.exists(disk):
        try:
            import json
            with open(disk) as fh:
                bases = set(json.load(fh))
            live = {canon.key(m) for m in host.metabolites
                    if _base(m.id).lower() in bases}
            if live:
                _PRODUCIBLE_CACHE[host] = (fp, live)
                report["status"] = "ok"
                return live
        except Exception:  # noqa: BLE001
            pass
    try:
        from cobra.flux_analysis import flux_variability_analysis
        # processes=1 is deliberate: cobra's multiprocessing FVA uses Windows spawn,
        # which re-executes the caller's __main__ in every worker — an unguarded script
        # (or a heredoc) fork-bombs. The disk cache makes this one-off cost acceptable.
        fva = flux_variability_analysis(host, fraction_of_optimum=0.0, processes=1)
    except Exception:  # noqa: BLE001
        report["status"] = "fva_failed"     # observable: the caller can warn the user
        return None
    live: Set[str] = set()
    live_bases: Set[str] = set()
    for m in host.metabolites:
        for r in m.reactions:
            c = r.metabolites[m]
            try:
                lo, hi = fva.at[r.id, "minimum"], fva.at[r.id, "maximum"]
            except Exception:  # noqa: BLE001
                lo, hi = r.lower_bound, r.upper_bound
            if (c > 0 and hi > 1e-9) or (c < 0 and lo < -1e-9):
                live.add(canon.key(m))
                live_bases.add(_base(m.id).lower())
                break
    if not live:
        report["status"] = "no_evidence"
        return None
    try:
        _PRODUCIBLE_CACHE[host] = (fp, live)
    except TypeError:
        pass
    try:
        import json
        with open(disk, "w") as fh:
            json.dump(sorted(live_bases), fh)
    except Exception:  # noqa: BLE001
        pass
    report["status"] = "ok"
    return live


def _native_keys(canon: "Canonicalizer", host: cobra.Model,
                 require_producible: bool) -> Set[str]:
    """The host compounds a route may bottom out in."""
    present = {canon.key(m) for m in host.metabolites}
    if not require_producible:
        return present
    live = producible_keys(host, canon)
    return present if live is None else (present & live)


def _name_tokens(s: str) -> Set[str]:
    import re as _re
    return {t for t in _re.split(r"[^a-z0-9]+", (s or "").lower()) if len(t) > 2}


def _ngrams(tokens: Set[str], n: int = 7) -> Set[str]:
    """All ``n``-character substrings of the given name tokens (chemical-root matching)."""
    out: Set[str] = set()
    for t in tokens:
        for i in range(0, max(0, len(t) - n + 1)):
            out.add(t[i:i + n])
    return out


def _inchikey_block_of(smiles: str) -> str:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(m).split("-", 1)[0] if m is not None else ""
    except Exception:  # noqa: BLE001
        return ""


def nearest_reachable_analogues(db: cobra.Model, target_name: str, *,
                                target_smiles: str = "", limit: int = 8) -> List[dict]:
    """Compounds in ``db`` that RESEMBLE the target and can actually be produced (L6).

    When a target cannot be found or reached, "not found" is the least useful possible
    answer. Very often a close relative *is* present and reachable — violacein is stored
    as ``violaceinate``; N,N-dimethyltryptamine is absent but ``N-methyltryptamine`` is
    there; THCV is absent but ``cannabidiolic acid`` is. Surfacing those turns a dead end
    into "here is what you CAN reach, and your target is one step further".

    Candidates are ranked by structure first (InChIKey connectivity block == the same
    skeleton), then by name-token overlap, and are required to have at least one
    producing reaction. Returns dicts with id/name/n_producers/score/reason.
    """
    want_block = _inchikey_block_of(target_smiles) if target_smiles else ""
    want_tokens = _name_tokens(target_name)
    # Chemical names share roots rather than whole words: "tryptamine" sits inside
    # "dimethyltryptamine", "cannabi" links "tetrahydrocannabivarin" to "cannabidiolic
    # acid", and "violacein" to "protoviolaceinate". Comparing sets of 7-character
    # n-grams catches all three, while staying specific enough to reject unrelated
    # neighbours ("butadiene" vs "butanoate" share only 4 characters).
    want_grams = _ngrams(want_tokens)

    # Not all shared roots are equally meaningful. "tetrahydro", "hydroxy" and "methyl"
    # decorate half the database and say nothing about identity, while "cannabi" or
    # "violacei" nearly pin the compound down. Weight each shared n-gram by how RARE it
    # is across the database (inverse document frequency), so a single distinctive root
    # outranks ten generic ones.
    import math
    gram_df: Dict[str, int] = {}
    candidates: List[tuple] = []
    for m in db.metabolites:
        producers = 0
        for r in m.reactions:
            c = r.metabolites[m]
            if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
                producers += 1
        if producers <= 0:
            continue                      # an orphan is not a useful suggestion
        toks = _name_tokens(m.name or "")
        grams = _ngrams(toks)
        candidates.append((m, producers, toks, grams))
        for g in grams:
            gram_df[g] = gram_df.get(g, 0) + 1
    n_docs = max(1, len(candidates))

    scored: List[dict] = []
    for m, producers, toks, grams in candidates:
        name = m.name or ""
        score, reason = 0.0, ""
        ann = m.annotation or {}
        blk = ann.get("inchi_key")
        if isinstance(blk, list):
            blk = blk[0] if blk else None
        if want_block and blk and str(blk).split("-", 1)[0] == want_block:
            score, reason = 1.0, "same molecular skeleton (InChIKey match)"
        elif want_tokens and toks & want_tokens:
            score, reason = 0.85, "shares a name word with the target"
        else:
            shared = want_grams & grams
            if shared:
                # Rank by the RAREST shared root: one distinctive root ("cannabi") is
                # far stronger evidence than many generic ones ("tetrahyd", "hydroxy").
                best = min(shared, key=lambda g: gram_df.get(g, n_docs))
                df = gram_df.get(best, n_docs)
                idf = math.log(n_docs / max(1, df)) / math.log(n_docs)   # 0..1
                if idf >= 0.45:            # discard ubiquitous decorator roots
                    score = min(0.8, 0.35 + 0.5 * idf)
                    reason = f"distinctive shared root “{best}” ({df} other compounds)"
        if score <= 0:
            continue
        scored.append({"id": m.id, "name": name or m.id, "n_producers": producers,
                       "score": round(score, 3), "reason": reason})
    scored.sort(key=lambda d: (-d["score"], -d["n_producers"], d["id"]))
    return scored[:limit]


def route_producibility_report(host: cobra.Model, db: cobra.Model,
                               reaction_ids: List[str]) -> dict:
    """Second-pass producibility check for an already-found route (never a search filter).

    The search deliberately does NOT gate on producibility — a route grounding on a
    host compound that carries zero flux under the current medium is still a valid
    *candidate* the user may want to enable. This is the on-demand FVA the user runs to
    check whether the route will actually flow, and it must WARN rather than hide.

    Returns ``{"status": ok|fva_failed|no_evidence|no_reactions,
    "nonproducible": [(host_id, name), …]}``. ``status != "ok"`` means the check could
    not be computed and the caller should say so, not silently assume the route is fine.
    """
    canon = Canonicalizer([host, db])
    currency = _currency_keys(canon, [host, db])
    present = {canon.key(m): m for m in host.metabolites}
    report: dict = {}
    live = producible_keys(host, canon, report)
    out = {"status": report.get("status", "ok"), "nonproducible": []}
    if live is None:
        return out                          # status explains why we can't judge
    seen: Set[str] = set()
    for rid in reaction_ids:
        if not db.reactions.has_id(rid):
            continue
        for m, c in db.reactions.get_by_id(rid).metabolites.items():
            if c >= 0:
                continue                    # only consumed precursors ground the route
            k = canon.key(m)
            if k in present and k not in live and k not in currency and k not in seen:
                seen.add(k)
                hm = present[k]
                out["nonproducible"].append((hm.id, (getattr(hm, "name", "") or hm.id)))
    return out


def _start_keys(canon: "Canonicalizer", host: cobra.Model, db: cobra.Model,
                start_metabolites: Optional[List[str]]) -> Optional[Set[str]]:
    """Canonical keys for a user-supplied starting-compound list, or ``None``."""
    if not start_metabolites:
        return None
    keys = set()
    for sid in start_metabolites:
        met = _resolve_target(host, db, sid)
        if met is not None:
            keys.add(canon.key(met))
    return keys or None


def retro_search(host: cobra.Model, db: cobra.Model, target_id: str, *,
                 max_reactions: int = 30, preferred_ec: Optional[set] = None,
                 forbidden: Optional[set] = None, include_boundary: bool = False,
                 require_balanced: bool = False,
                 start_metabolites: Optional[List[str]] = None,
                 require_producible: bool = False) -> Optional[List[str]]:
    """Minimal set of heterologous reactions connecting ``target_id`` back to
    compounds the host already makes. Precursor-first order. ``None`` if unreachable.

    ``start_metabolites`` *restricts* the set of precursors a route may bottom out in
    (default: everything the host makes). It is a narrowing constraint and must never
    widen what is available. ``require_producible`` additionally demands that those
    precursors be compounds the host can actually produce, not merely ones present in
    the model (see :func:`producible_keys`).
    """
    from .databases import reaction_ec_numbers
    from .pathway_design import _is_grossly_unbalanced

    forbidden = forbidden or set()
    preferred_ec = preferred_ec or set()
    canon = Canonicalizer([host, db])
    currency = _currency_keys(canon, [host, db])
    starts = _start_keys(canon, host, db, start_metabolites)
    base_native = (starts if starts is not None
                   else _native_keys(canon, host, require_producible))
    native = base_native | currency

    target_met = _resolve_target(host, db, target_id)
    if target_met is None:
        return None
    target_key = canon.key(target_met)
    if target_key in native:
        return []          # host already makes it — no heterologous reactions

    # producers[compound_key] = list of (pref, rid, substrate_keys)
    producers: Dict[str, list] = {}
    # Reactions that verifiably do not balance. They are NOT dropped — a missing
    # formula or an odd protonation state must not hide real chemistry — but a route
    # through one must never be preferred over an equally short balanced route.
    unbalanced: Set[str] = set()
    for r in db.reactions:
        if r.id in forbidden:
            continue
        if not include_boundary and r.boundary:
            continue
        gross = _is_grossly_unbalanced(r)
        if require_balanced and gross:
            continue
        if gross:
            unbalanced.add(r.id)
        pref = 1 if (preferred_ec and set(reaction_ec_numbers(r)) & preferred_ec) else 0
        subs = {canon.key(m) for m, c in r.metabolites.items() if c < 0}
        prods = {canon.key(m) for m, c in r.metabolites.items() if c > 0}
        sub_need = {s for s in subs if s not in currency}
        prod_need = {p for p in prods if p not in currency}
        # A compound on BOTH sides is not *made* by this reaction. Skipping those is
        # what stops a compound being conjured from nothing: this search is
        # compartment-agnostic, so a pure transport (A_c <=> A_e) collapses to A <=> A.
        # Subtracting the product from its own substrate list then left an EMPTY
        # requirement, making every compound that has a transport reaction anywhere in
        # the database freely available in one step — which is how routes came to rely
        # on an "external supply" of non-native compounds such as hexanoyl-CoA. The
        # same guard correctly excludes catalytic/recycled carriers.
        for p in prod_need - sub_need:
            producers.setdefault(p, []).append((pref, r.id, sub_need))
        if r.reversibility or (r.lower_bound < 0 < r.upper_bound):
            for s in sub_need - prod_need:
                producers.setdefault(s, []).append((pref, r.id, prod_need))
    # Deterministic candidate order. Preferred-EC reactions first, then a stable
    # tie-break on the reaction id: `subs`/`prods` are sets of strings, whose iteration
    # order changes between processes (hash randomisation), so without this the same
    # query could return a different route — or no route at all — on each run.
    for lst in producers.values():
        lst.sort(key=lambda t: (-t[0], t[1]))

    # Bottom-up fixed point over the AND/OR graph: repeatedly relax every producer
    # until no route improves. A compound only acquires a route once ALL its
    # substrates already have one, so routes are precursor-first by construction and
    # cycles simply never close — no recursion, no cycle cut, no memo.
    #
    # This replaced a recursive depth-first `produce()` with a `memo[c] = None` cache.
    # That cache was unsound: a `None` from hitting a compound already on the current
    # branch is a failure only for THAT branch, but it was cached globally and poisoned
    # every later branch. Combined with iterating sets of strings (whose order changes
    # per process under hash randomisation), the same query could return a route, a
    # DIFFERENT route, or no route at all on successive runs. Relaxation to a fixed
    # point is order-independent, so the answer is now stable and reproducible.
    routes: Dict[str, List[str]] = {}
    best_key: Dict[str, tuple] = {}

    def route_for(c: str) -> Optional[List[str]]:
        return [] if c in native else routes.get(c)

    for _round in range(max(2, max_reactions) + 2):
        changed = False
        for c, plist in producers.items():          # dict order == db order: stable
            if c in native:
                continue
            for pref, rid, subs in plist:
                acc: List[str] = []
                ok = True
                for s in sorted(subs):
                    sr = route_for(s)
                    if sr is None:
                        ok = False
                        break
                    for x in sr:
                        if x not in acc:
                            acc.append(x)
                if not ok:
                    continue
                if rid not in acc:
                    acc.append(rid)
                if len(acc) > max_reactions:
                    continue
                # Fewest steps wins; then FEWEST verifiably unbalanced steps; then
                # preferred-EC; then route content for determinism.
                #
                # The balance term matters more than it looks. Merging databases brought
                # in `R12488` ("H2O + 4-coumaroyl-CoA <=> 4-hydroxyphenyllactate"),
                # which loses the entire CoA moiety — 21 carbons. It offered a 3-step
                # route to resveratrol, exactly as long as the real tyrosine
                # ammonia-lyase route, so the tie was settled by the LAST term:
                # "34HPLFM" sorts before "L_tyrosine_ammonia_lyase", and the engine
                # returned chemically impossible chemistry (yield 3.6e-05 instead of
                # 0.29) purely because of alphabetical order.
                n_bad = sum(1 for x in acc if x in unbalanced)
                key = (len(acc), n_bad, -pref, tuple(acc))
                if c not in best_key or key < best_key[c]:
                    routes[c], best_key[c] = acc, key
                    changed = True
        if not changed:
            break

    return route_for(target_key)


# --- forward expansion search ------------------------------------------------
def expansion_search(host: cobra.Model, db: cobra.Model, target_id: str, *,
                     max_steps: int = 20, preferred_ec: Optional[set] = None,
                     forbidden: Optional[set] = None, include_boundary: bool = False,
                     require_balanced: bool = False,
                     start_keys: Optional[set] = None,
                     start_metabolites: Optional[List[str]] = None,
                     require_producible: bool = False) -> Optional[List[str]]:
    """Forward network expansion over the identity space. Fires reactions whose
    non-currency substrates are all reachable until the target is produced.

    ``start_metabolites`` (ids, resolved here) or ``start_keys`` (canonical keys)
    *restrict* the starting pool; both are narrowing constraints.
    """
    from .databases import reaction_ec_numbers
    from .pathway_design import _is_grossly_unbalanced

    forbidden = forbidden or set()
    preferred_ec = preferred_ec or set()
    canon = Canonicalizer([host, db])
    currency = _currency_keys(canon, [host, db])
    starts = set(start_keys) if start_keys else _start_keys(canon, host, db, start_metabolites)
    reachable = (set(starts) if starts
                 else _native_keys(canon, host, require_producible))
    reachable |= currency

    target_met = _resolve_target(host, db, target_id)
    if target_met is None:
        return None
    target_key = canon.key(target_met)
    if target_key in reachable:
        return []

    directions = []  # (pref, bad, rid, sub_keys, prod_keys)
    for r in db.reactions:
        if r.id in forbidden:
            continue
        if not include_boundary and r.boundary:
            continue
        gross = _is_grossly_unbalanced(r)
        if require_balanced and gross:
            continue
        pref = 1 if (preferred_ec and set(reaction_ec_numbers(r)) & preferred_ec) else 0
        subs = {canon.key(m) for m, c in r.metabolites.items() if c < 0}
        prods = {canon.key(m) for m, c in r.metabolites.items() if c > 0}
        sub_nc = {s for s in subs if s not in currency}
        bad = 1 if gross else 0
        directions.append((pref, bad, r.id, sub_nc, prods))
        if r.reversibility or (r.lower_bound < 0 < r.upper_bound):
            prod_nc = {p for p in prods if p not in currency}
            directions.append((pref, bad, r.id, prod_nc, subs))
    # Preferred-EC first, then BALANCED reactions ahead of verifiably unbalanced ones
    # (they are kept, but must never win by luck — see the note in retro_search), then
    # a stable tie-break on the reaction id so the firing order is reproducible.
    directions.sort(key=lambda d: (-d[0], d[1], d[2]))

    pred: Dict[str, tuple] = {}
    found = False
    for _ in range(max(2, max_steps) + 6):
        changed = False
        for _pref, _bad, rid, sub_nc, prods in directions:
            if sub_nc <= reachable:
                for p in sorted(prods):        # sorted: deterministic across runs
                    if p in currency or p in reachable:
                        continue
                    reachable.add(p)
                    pred[p] = (rid, sub_nc)
                    changed = True
        if target_key in reachable:
            found = True
            break
        if not changed:
            break
    if not found:
        return None

    ordered: List[str] = []
    seen = set()
    stack = [target_key]
    visited = set()
    while stack:
        k = stack.pop()
        if k not in pred or k in visited:
            continue
        visited.add(k)
        rid, subs = pred[k]
        if rid not in seen:
            seen.add(rid)
            ordered.append(rid)
        for s in sorted(subs, reverse=True):   # sorted: deterministic across runs
            if s not in visited:
                stack.append(s)
    ordered.reverse()
    return ordered


# --- database merging --------------------------------------------------------
def merge_databases(models: List[cobra.Model], *, name: str = "merged_universal"):
    """Unify several databases into one cytosolic model with each compound and
    reaction present exactly once (dedup by identity / stoichiometry)."""
    from .pathway_design import readable_metabolite_id, readable_reaction_id

    models = [m for m in models if m is not None]
    if not models:
        return cobra.Model(name)
    canon = Canonicalizer(models)
    merged = cobra.Model(name)
    met_by_key: Dict[str, cobra.Metabolite] = {}
    used_ids: Set[str] = set()

    def merged_met(met: cobra.Metabolite) -> cobra.Metabolite:
        k = canon.key(met)
        if k in met_by_key:
            m = met_by_key[k]
            # enrich annotation from duplicates so future matches are richer
            for a, v in (getattr(met, "annotation", {}) or {}).items():
                m.annotation.setdefault(a, v)
            return m
        readable = readable_metabolite_id(met)
        base = _base(readable)
        cid = f"{base}_c"
        i = 2
        while cid in used_ids:
            cid = f"{base}_{i}_c"
            i += 1
        used_ids.add(cid)
        nm = cobra.Metabolite(cid, name=met.name or base, formula=met.formula,
                              charge=met.charge, compartment="c")
        nm.annotation = dict(getattr(met, "annotation", {}) or {})
        met_by_key[k] = nm
        return nm

    seen_sig: Set[tuple] = set()
    seen_rxn_ids: Set[str] = set()
    new_reactions = []
    for model in models:
        for rxn in model.reactions:
            if rxn.boundary:
                continue
            coeffs: Dict[cobra.Metabolite, float] = {}
            for met, c in rxn.metabolites.items():
                mm = merged_met(met)
                coeffs[mm] = coeffs.get(mm, 0.0) + c
            coeffs = {m: c for m, c in coeffs.items() if abs(c) > 1e-9}
            if len(coeffs) < 2:                     # transport/identity after collapse
                continue
            sig = tuple(sorted((m.id, round(c, 6)) for m, c in coeffs.items()))
            if sig in seen_sig:
                continue
            seen_sig.add(sig)
            rid = readable_reaction_id(rxn)
            base_rid, j = rid, 2
            while rid in seen_rxn_ids:
                rid = f"{base_rid}_{j}"
                j += 1
            seen_rxn_ids.add(rid)
            nr = cobra.Reaction(rid, name=rxn.name, subsystem=getattr(rxn, "subsystem", ""))
            nr.bounds = rxn.bounds
            nr.annotation = dict(getattr(rxn, "annotation", {}) or {})
            nr._pending = coeffs
            new_reactions.append(nr)

    merged.add_metabolites(list(met_by_key.values()))
    merged.add_reactions(new_reactions)
    for nr in new_reactions:
        nr.add_metabolites(nr._pending)
        del nr._pending
    return merged
