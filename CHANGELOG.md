# Changelog

All notable changes to **GSM ToolBox** are recorded here.

## 0.3.10 - 2026-08-10

### Added: reconciling metabolite identifiers across a model and a database

A model and a reaction database can use different identifiers for the same compound.
Where their cross-references overlap the toolbox already unified them. Where the
annotation sets are disjoint it could not, and the database metabolite entered the model
as an orphan. The route then carried no flux while looking complete, with no error shown.

`core/id_reconcile.py` proposes correspondences from several kinds of evidence, each
scored and each reported: InChIKey, shared cross-references, formula and charge,
normalised names. Names must match in full after normalisation, because Zeaxanthin is a
substring of Zeaxanthin diglucoside and those are different compounds. Nothing is applied
automatically below a high threshold, and neither model is modified.

Protein-bound cofactors are treated separately through a curated allow-list. Ferredoxin
is written `Fe2S2X` by one database and `Fe2R8S2` by another, and neither is a real
molecular formula, so comparing them elementally rejects the very match that is needed.
For those classes the formula is recorded as carrying no evidence. Redox state is still
enforced, so reduced ferredoxin can never merge into the oxidised form.

*Tools > Pathway Design > Reconcile metabolite identifiers* lists the proposals with
their score and evidence. Accepted mappings are stored with the project, versioned, and
listed in any report the project produces. A correspondence the checks rejected can still
be applied, and is recorded as an override.

### Fixed: a route to a real target carried no flux

Routing fucoxanthin through a cyanobacterial model needed a database reaction consuming
reduced ferredoxin. The generic cofactor resolver missed it twice over: the database
writes the name in the plural where the alias table is singular, and the host spells the
oxidised form `fdxo_2_2`, which no alias listed. The lookup is now plural insensitive
and
falls back to matching the cofactor class, so a host that spells a carrier differently no
longer needs a new entry added by hand. The route now carries flux through the model's
own native carotenoid backbone rather than re-adding it.

### Changed: documentation

Installation instructions now cover a portable Windows build alongside the installer, and
the macOS and Linux routes are described in full. Author and licence information credits
Jorge Fernandez Mendez throughout. The manual and the shipped documents were rewritten in
plainer prose.

## 0.3.9 - 2026-08-17

### Fixed. A metabolite with no formula was reported as a mass imbalance
`Fucoxanthin_c` has no formula in the merged database, so COBRApy's
`check_mass_balance()` counted its atoms as **zero** and reported
`Fucoxanthin_synthase` as `{'C': -42, 'H': -58, 'O': -6}`. The reaction accused of
losing 42 carbons when nothing is wrong with it and the product's atoms simply were
never counted. Fetching the formula (`C42H58O6`) shows the reaction to be perfectly
**balanced**. Roughly **40% of the merged database's 26,322 metabolites carry no
formula**, so this was a systematic false accusation, not an edge case.

Every balance check now has **three** outcomes rather than two. Balanced, imbalanced,
and *cannot be checked*. Carried by a single `BalanceVerdict` type and reported
separately everywhere: the reaction editor and builder, the pathway result and its
warnings, the feasibility report, the apply-pathway summary, the route comparison, and
the model quality report. The uncheckable case **names the participant to blame**
("cannot be checked. Fucoxanthin (Fucoxanthin_c) has no formula") instead of leaving
a bare "unknown". The quality score no longer docks a model for formulas the database
never shipped.

### Fixed. Carbon yield reported `nan` instead of saying why
A route whose product has no formula produced a bare NaN, which reads as a failed
design rather than a fixable data gap. Carbon yield now reports **"not computable"**
with the reason: the product has no formula, it contains no carbon, the route carries
no flux, or no carbon uptake was resolved. A *consumed* metabolite with no formula now
withholds the figure too, rather than under-counting the denominator and inflating the
yield. The same absent-versus-zero error hidden on the other side of the fraction.

### Added. Fetch a missing formula on demand
Where a check cannot run, the app now offers to close the gap instead of describing
it. **Fetch missing formula** in the reaction editor, and on the right-click menu of an
uncheckable pathway step, resolves the formula from BiGG's own table first (the same
source as the stoichiometry), then the BiGG API, then the compound's
InChIKey/KEGG/ChEBI cross-references, and re-runs the check. A formula that cannot be
resolved is left absent, never guessed, which would turn an honest "cannot be checked"
into a confident wrong answer. An absent charge likewise stays absent.

## 0.3.8 - 2026-08-09

### Fixed. Series colours on dynamic plots were not distinguishable
The plot palette opened with two blues, `#2D6CDF` and `#1A73E8`, **2.1 apart in
OKLab**. One colour as far as a line plot is concerned, and held only seven entries,
so an eighth reaction repeated the first. Following several fluxes at once produced
lines that could not be told apart.

The palette is now **twelve colours chosen by maximising the worst pairwise
separation**, not by taste: worst adjacent pair ΔE 28 for full-colour vision and 11.3
under simulated protanopia. Twelve mutually distinguishable colours is achievable for
normal vision; twelve that also survive red-green colour blindness is not, because
protanopia and deuteranopia collapse the red-green axis. Plots therefore also vary
**marker shape** per series and **dash pattern** past the twelfth, so identity never
rests on hue alone. The distances are asserted in the test suite rather than claimed.

Also fixed: with many series the y-axis label concatenated every column name and ran
off the figure. Past three series the axis defers to the legend.

### Added. Plot Results inside the run's table window
The results popup for a stored scan or time course now carries a **Plot Results**
button, so the plot is reachable from where the numbers are being read.

## 0.3.7 - 2026-08-09

### Fixed. Escher Visualizer opened nothing
`EscherExplorer` was constructed but never added to the tab widget, dropped from
`addTab` when the Dynamic Analysis tab was introduced. `setCurrentWidget` is a silent
no-op on a widget the tab bar does not own, so the map was built. Several seconds of
frozen interface. Into something that could never be shown. It is now inserted on
demand, with a wait cursor while the map builds. A test guards the whole class of
defect: a panel that exists but is unreachable.

### Fixed. A merged database could vanish while still on disk
The database list was built only from `db_registry.json`. When that file is missing, never written, deleted, or lost when the folder moves between machines. A merge that
cost hours disappeared from the interface while still occupying its space on disk.
Discovery now scans the databases folder directly, so the registry is a convenience
rather than the source of truth. Build intermediates, chemistry lookups and bundled
example models are filtered out; small fetched databases are not (a single fetched
pathway can be well under a kilobyte and is still a database the user asked for).

### Fixed. Green now means one thing on every network map
Green is reserved model-wide for metabolites that are **not native** to the model. The
peroxisome compartment was itself green and collided with it; hashed hues for unfamiliar
compartments could land in the green band by chance. Both paths now rotate clear of it,
rotating rather than clamping so two compartments cannot collapse onto one colour. A
legend above the map states what every colour means **for the model currently open**,
listing only compartments that actually carry metabolites.

### Changed. Dynamic Analysis rebuilt around stored runs
- The settings take the whole panel. The introductory prose and the regulation controls
  that sat above them are gone.
- One action row: **Run**, **Plot Results**, **Display Table**, a **Regulation**
  checkbox and its settings button. The Run button names the analysis it will start.
- Completed runs accumulate as closable, renamable tabs (`ConditionScan_1`,
  `TimeCourse_1`), each opening in its own non-modal window. Runs no longer overwrite
  one another. The comparison between two runs is usually the reason for the second.

### Added. Overlay plots across several runs
The plot dialog takes any number of stored runs and draws them on shared axes. Colour
is the quantity and line style is the run. Axis options are the union across the
selected runs, and a run missing a column is reported as missing rather than drawn as a
line at zero. Export is long form, one row per run per point.

### Added. Warning when a step swaps a compound for its isomer
ModelSEED `rxn35740`, "alcohol dehydrogenase [NAD(P)+]", pairs **n-butanol** with
**isobutanal**: a redox *and* a rearrangement of the carbon skeleton in one step, which
no alcohol dehydrogenase performs. It is mass- and charge-balanced, so every numeric
check passed it, and neither compound carries a structure, so the existing structural
skeleton check returned "unresolved" and said nothing.

Skeleton mismatches are now also detected from **nomenclature**. `n-` is straight,
`iso`/`2-methyl` is branched, which needs no structure and no network. Such a finding
is explicitly a suspicion: it reads "Possible", carries a `from_names` flag, and never
silently drops a route. True isomerases (identical formulas) stay silent.

### Changed. User-facing text
Illustrative examples naming particular compounds, organisms or genes have been removed
throughout; the explanatory text is shorter. Format placeholders that show what a field
expects (reaction-equation and gene-rule syntax) are deliberately kept. They document a
format rather than a case. A test keeps both halves of that distinction honest.

## 0.3.3 - 2026-07-22

### Added. Selenzyme-style enzyme search by reaction similarity
Finds candidate enzymes for a reaction **even when it has no EC number**. The case for
every rule-generated (RetroRules) step, and the capability 0.3.2 explicitly lacked.

- Enable once in **Settings ▸ Preferences ▸ Enable reaction-similarity enzyme search**.
  A one-time ~20 MB download of open **Rhea** data builds a 13 MB local index in ~90 s;
  everything then runs offline. Then use *Find enzymes by reaction similarity…* in the
  **Suggest EC numbers & enzymes** dialog.
- **8,711** characterised Rhea reactions indexed, joined to **397,133** UniProt enzyme
  mappings. Hits are ranked by reaction similarity, then curated (SwissProt) status and
  taxonomic proximity to your host; the closest-match similarity is reported so a loose
  match (< 0.35) is visibly flagged.

**Why it is a native implementation.** Upstream Selenzyme could not be installed: it
ships **only as a Docker service** with **2017-era pinned dependencies** (RDKit 2017,
Python 3.6) that cannot coexist with this app's RDKit 2026 / Python 3.11; there is **no
PyPI package**; its algorithm repository (`pablocarb/selenzy`) declares **no licence**, so
its code cannot be vendored; and the hosted service **refused connections**. The method is
published and its inputs are open data, so it is implemented directly on Rhea + UniProt.
Full rationale in `docs/enzyme_selection.md`.

**Two fixes that decided whether it worked** (both found by validation, not assumed):
- **Difference fingerprints, not structural ones.** A structural reaction fingerprint is
  dominated by molecular bulk and could not discriminate. The first build returned
  *ammonium transporters* for ethanol oxidation, and a nonsense decoy scored 0.172. The
  difference fingerprint drops that decoy to **0.010**.
- **Cofactors are stripped from both sides.** Written in full the reaction is mostly NAD,
  so a cofactor-free query (exactly how RetroRules emits steps) scored **0.089** against
  the real alcohol-dehydrogenase reaction. Stripping ubiquitous cofactors from reference
  and query leaves the transformation, and the same query then scores **1.000**.

**Validated:** `CCO>>CC=O` → **1.000, Alcohol dehydrogenase (EC 1.1.1.1)**;
tryptamine → N-methyltryptamine (the study's DMT step, which has no EC) → **0.674,
tryptophan methyltransferase *trpM* / Psilocybin synthase (EC 2.1.1.-)**; nonsense
chemistry correctly rejected.

### Fixed
- Rhea EC lookups now fall back from a directional reaction id to its master id
  (`rhea-reaction-smiles` is directional, `rhea2ec` is keyed on the master), so the EC
  column is populated instead of always blank.

## 0.3.2 - 2026-07-22

### Fixed. The "Namespace match. Please review" popup (a real mis-identification)
- A warning naming ambiguous matches such as **`CO2_c → co2_c, cobalt2_c`** appeared on
  almost every pathway. The cause was a **contaminated database entry**, not a display
  bug: the merged universal's `CO2` carries cobalt(2+) identifiers (KEGG `C00175`, ChEBI
  `48827`/`48828`/`23337`, HMDB `HMDB00608`, SEED `cpd00149`, InChIKey
  `XLJKHNWPARRRJB`). BioCyc writes cobalt as `CO+2`, which a merge fused with `CO2`.
- Worse, matching stopped at the **first cross-reference from an unordered set**, so
  carbon dioxide could silently be mapped onto **cobalt**.
- Matching now ranks candidates by **chemical identity first**. Elemental formula
  (decisive here: `CO2` vs `Co`), then InChIKey skeleton, then the number of supporting
  cross-references, and is deterministic. The warning now fires **only on a genuine
  tie**. Across the whole 26 322-metabolite merged database the ambiguous count went
  from constant to **zero**, and `CO2_c` maps correctly to `co2_c`.

### Added. EC numbers for reactions that don't carry one
- Right-click any reaction in a predicted pathway ▸ **"Suggest EC numbers & enzymes…"**.
  EC numbers are gathered from three kinds of evidence and each is labelled with its
  source: **annotated**, **cross-reference** (the reaction's own KEGG/Rhea id is looked
  up live), or **inferred** (a database reaction with an identical participant set
  carries one. Flagged *verify*). Every suggestion links to **ExPASy ENZYME** and
  **BRENDA**. Example: `R05380` has no EC annotation but resolves to **EC 4.2.1.112**
  via KEGG.

### Added. Enzyme candidates ("clone this gene")
- From a suggested EC, **Find enzyme candidates (UniProt)** lists sequences ranked with
  **curated (SwissProt) entries first** and organisms **taxonomically related to your
  production host** promoted (cyanobacteria, enterobacteria, yeasts, *Pseudomonas*,
  *Bacillus*, actinobacteria). Ranking never hides a candidate. Double-click opens the
  UniProt entry.
- **Selenzyme:** it was never implemented. It existed only as a plan in
  `pathway_design_strategy.md` §2.3. Its public service was probed and **refused the
  connection**, and it has no versioned public API, so the capability is implemented on
  UniProt/KEGG/Rhea instead of depending on it. Full rationale and limitations in
  **`docs/enzyme_selection.md`**.

### Added. EC links in RetroRules
- Rule-based steps render their EC numbers as **clickable ExPASy links**, so a proposed
  step leads straight to the enzyme class that performs it.

## 0.3.1 - 2026-07-22

Findings from an eleven-compound production study (see
`docs/production_study/`) drove this release: the study exposed a false-negative in target
resolution, a reproducibility hole in rule search, and several numbers that were being
reported more confidently than they deserved.

### Changed. Thermodynamics (MDF) is now opt-in
- The MDF suite is **hidden by default** and enabled from **Settings ▸ Preferences ▸
  Enable MDF Suite**. Ticking it asks for consent and downloads the eQuilibrator
  compound database (~1.34 GB) once; until then no MDF control appears anywhere.
- **The compound database is no longer bundled**, superseding 0.3.0. It tripled the
  installer for a feature most users never touch and presented an assumption-laden
  analysis as a headline capability. Only the (small) eQuilibrator code ships.
- MDF results now state the physiology they rest on (pH 7.5, I 0.25 M, pMg 3.0,
  1 µM–10 mM, 298 K), and a **single-reaction route is flagged**: its driving force is
  not a pathway MDF and its magnitude is not comparable with multi-step routes.

### Fixed. Reproducibility and honest reporting
- **Rule search is deterministic (L1).** Rules are explored in a stable, sorted order
  permuted by a user-visible **seed** (Preferences). The same query now returns the same
  routes every run. Verified on butane, which previously flipped between a 3.46 flux and
  zero. The RetroRules dialog states the seed and that results are a *sample*, not an
  exhaustive enumeration.
- **Rule-route flux is no longer reported as a capacity (L2).** It was measured on a
  synthetic demand for a novel compound and saturated against a generic host bound. Four unrelated products returned the identical 3.46. Such routes now report only
  whether they carry flux.
- **Branching yields are no longer null (L11).** The product sink is resolved across
  compartments (following the transport to the *_e* exchange), the same gap that made
  capacity calculations silently return zero for exported products.

### Added
- **Carbon yield (L10)** alongside raw flux. Mol C in product / mol C consumed. The
  figure that IS comparable between targets of different size (ethylene: 72.7%).
- **Find strategies (FSEOF)** button on any route, running on the host-integrated model,
  so **rule-based routes get strategy analysis too (L9)** rather than only database ones.
- **Nearest reachable analogue (L6).** When no route is found, the app now suggests
  related compounds that *can* be produced, ranked by shared distinctive chemical roots
  (rarity-weighted) and requiring a producing reaction. This is what turns "violacein is
  not in the database" into "violaceinate is, with 9 producing reactions", and surfaces
  N-methyltryptamine for DMT and cannabidiolic acid for THCV.
- **Add a final step (L6).** Declare the one known enzyme between a reachable compound
  and the real target, and the completed route is built and evaluated end-to-end.

### Fixed. RetroRules bookkeeping
- **Versioned rule cache (L12).** The per-diameter pickle now carries a schema version,
  so adding a parsed field (as `Score_normalized` was) forces a rebuild instead of
  silently serving rules with score 0 and degrading "rank by reliability" to step count.
  Stale caches are purged.
- **Readable identifiers (L13).** Rule reactions are now `RR_tryptamine_synthesis`
  rather than `RR_step2`, so they are traceable in a model, table or report.

## 0.3.0 - 2026-07-22

### Added. Thermodynamic feasibility (MDF)
- A **Thermodynamics (MDF)** button appears on any found pathway. It computes each
  reaction's ΔrG′° via **eQuilibrator** and solves the pathway **Max-min Driving Force**
  (Noor et al. 2014): a positive MDF means concentrations exist that drive every step
  forward; otherwise the least-favourable step is named as the bottleneck. The MDF linear
  program is SciPy-only and unit-tested; eQuilibrator and its ~1.3 GB compound cache are
  bundled so this works offline (the cache is placed on first use, or downloaded on demand
  if a build ships without it).

### Added. Branching / competition analysis
- A **Branching & competition** button surfaces the EA-MNE analysis on any route: ideal
  (linear) vs realistic (network) yield, the percentage lost to competing reactions, and
  the competing reactions themselves as one-click knockdown candidates.

### Added. More chemistry sources
- **Fetch missing chemistry** now offers **Rhea** (open, ChEBI-based) alongside KEGG, and
  explains how to reach licence-restricted **MetaCyc** content (local BioCyc file, or via
  the MetaNetX universal).

### Added. Update checker
- **Help ▸ Check for updates** reports the installed version of the toolbox, its key
  dependencies, and its downloadable data (databases, RetroRules ruleset, eQuilibrator
  cache), flags anything out of date against PyPI, and suggests the manual command to
  update each. It never installs anything itself.

### Changed. Producibility is a warning, not a filter
- The FVA producibility check no longer prevents a route from being found. Every candidate
  route is shown; the check runs as a **second pass** (in *Flux feasibility*) that warns
  when a route grounds on a host compound that carries zero flux under the current medium.
  Its result is cached to disk (keyed by a model fingerprint) so the FVA is paid once, and
  it now runs single-process to avoid the Windows multiprocessing fork-bomb. A failed FVA
  is reported instead of silently trusting everything.

### Improved. RetroRules usability
- The RetroRules dialog now shows the **target's structure and name at the top**, offers
  **several ranked alternative routes** (ranked by RetroRules reliability score, fewest
  steps, or fewest native precursors. User's choice), and lets you **tick which steps to
  add** instead of all-or-nothing. A new options dialog sets the number of alternatives,
  depth, and ranking before searching.

### Improved. Escher FBA navigator
- A one-line hint bar now appears in navigator mode ("click any reaction on the map to
  edit its flux"), making the click-to-edit affordance discoverable.

### Fixed. Window width
- The new pathway-panel action buttons wrap in a `FlowLayout`, keeping the window within a
  1366-px screen (guarded by a regression test).

## 0.2.5 - 2026-07-21

### Improved. RetroRules naming
- Compound names in the suggestions dialog and the added COBRA metabolites are now
  resolved much more robustly: the loaded databases first (by InChIKey), then an
  **online PubChem lookup** (common synonym, else IUPAC name, casing normalised) for the
  rule-generated intermediates the databases don't contain, so far fewer compounds show
  a bare SMILES.
- Added reactions carry a readable name (**"<product> synthesis"**), and the **suggested
  reaction id now defaults to that name** (editable later in the Apply dialog).

### Fixed. Diagnostics
- "Why not found" now names the blocking compound for **circular precursor pairs**
  (e.g. cyclopentanol ⇌ cyclopentanone, which are isolated from central metabolism in
  every database) instead of leaving the reason blank.

## 0.2.4 - 2026-07-20

### Changed. Pathway Design layout
- The three database buttons (**Manage / Load / Fetch**) moved to the top of the
  right-hand *Reaction databases* panel, freeing the top row for a **wider target
  selector** beside **Predict pathway** and **RetroRules Prediction**. Trailing “…”
  removed from the button labels.

### Improved. RetroRules GUI
- The target's structure is now resolved automatically: SMILES/InChI annotation first,
  then an **online fetch** (PubChem by InChIKey/name, then KEGG) for compounds that carry
  only an id, no more asking for a SMILES except as a last resort.
- The suggestions dialog now **draws the 2-D structure** of every compound in each step
  (substrates → product) and labels them with **human-readable names** where the loaded
  databases recognise the compound.
- **Add selected reactions** turns a rule-based route into a suggested pathway in the
  main panel, so the normal **Apply pathway** flow (with the reaction/metabolite renaming
  popup) can add it to the model.

### Fixed. Diagnostics
- “Why not found” now **names the blocking compound for a circular pair**. Cyclopentanol
  ⇌ cyclopentanone each ‘produce’ the other but neither is reachable from central
  metabolism in any database (a genuine data gap. RetroRules is the route for it); the
  tool now reports *“depends on Cyclopentanone, which no reaction in this database
  produces”* instead of leaving the reason blank.

## 0.2.3 - 2026-07-20

### Added. RetroRules in the GUI
- A **Rule-based (RetroRules)** button in Pathway Design. It resolves the target's
  structure (SMILES/InChI, or asks for a SMILES), offers to download the full RR02
  dataset the first time (~43 MB), searches reaction rules behind a cancelable progress
  dialog, and shows the suggested precursor-first steps. Each labelled as a
  **prediction to verify**, distinct from database-backed reactions. The control row
  wraps, so the button cannot push the window past the screen.

## 0.2.2 - 2026-07-15

Interface fixes from testing, plus the pathway engine learning to **explain itself**
instead of reporting a bare `0`.

### Fixed. Interface
- **The window no longer leaves its maximized state.** The real cause was the status
  bar: its label reported the full width of whatever message it held, so a long
  completion message ("Saved strategy 'Round 2. …'", "FBA complete: …") pushed the
  window's minimum past the screen, and Qt honours a minimum size over the maximized
  state. Since almost every analysis posts a message, it fired constantly. Messages are
  now elided (full text on hover) and the status bar's minimum is constant.
- **Escher**: information moved to **right-click ▸ Show Information**, leaving
  left-click free to select a reaction and edit its flux; arrowheads are redrawn after
  a node is dragged, so they no longer keep a stale angle; the value label reads
  **Flux** rather than "Data".
- **Strategy Visualizer**: **Run FBA** / **Run pFBA** without leaving the tab, and the
  reaction to follow can now be **any reaction**. Previously only exchanges were
  offered, so a heterologous product (which often has no exchange) could not be plotted.
- **Pathway map**: metabolites taken from the database are drawn **green**, host-native
  precursors **blue**, so what is made de novo is obvious at a glance.
- A model whose external compartment cannot be guessed no longer breaks the Strategy
  tab (cobra's `.exchanges` raises on such models).

### Added. Pathway diagnostics
- **"Flux & yield analysis"** on any route. It separates the two causes of zero flux
  that look identical to the user: a route that is genuinely **blocked** (naming the
  bottleneck step and the dead-end compound), and one that is perfectly capable but
  **not rewarded by the objective**. The usual reason a design predicting 0.21 shows 0
  after FBA, which maximises biomass rather than the product.
- **Competing-reaction (branching) analysis**, after EA-MNE (*J. Chem. Inf. Model.*
  2026): reports which reactions divert the route's intermediates, what that costs in
  yield, and distinguishes branches that can be knocked out from essential ones that
  can only be down-regulated.
- **"Find a route that runs"** retries the search around the blocking step.
- **Zero-flux routes are offered, not hidden**. With the reason attached.

### Added. "why not found"
`0 pathways` is no longer an answer. The search now distinguishes: the compound is not
in the database; it is catalogued but **nothing produces it**; or a **named precursor**
is unreachable, and only the last is helped by raising the step limit, which is why
doing so so often changed nothing. Where a gap is named, **Fetch missing chemistry**
pulls KEGG reactions around *that compound*.

### Fixed. Engine
- **Generic placeholders understood.** ModelSEED's `NAD-P-OR-NOP`/`NADH-P-OR-NOP` gate
  **363 reactions each** and were treated as compounds needing synthesis, making every
  one of them unusable. Also generic amino-acid classes and protein-bound residues
  across all three databases.
- **A balanced route now beats an equally short unbalanced one.** Merging databases
  surfaced a reaction that hydrolyses 4-coumaroyl-CoA straight to
  4-hydroxyphenyllactate. Losing the entire CoA (21 carbons). It gave a 3-step route to
  resveratrol as short as the real one, and the alphabetical tie-break picked the
  impossible chemistry (yield 3.6e-05 vs 0.29). Ranking now prefers balanced routes on a
  tie; unbalanced reactions are still offered when they are the only option.
- Merging databases is much faster: the canonical key is memoised instead of being
  re-derived once per metabolite per reaction.
- The KEGG "focused database" fetch named the seed compound first with an adequate
  budget: previously the naming pass gave up after 15 s, so fetching "erythromycin"
  returned a database where the target was labelled `C01912` and unfindable by name.

### Added. Databases
- **Merging warns first, saves after.** A merge states its expected cost up front (up to
  hours for the big three), keeps your source databases instead of replacing them, and
  **saves + registers the result** so it loads instantly next time. Databases can be
  **renamed** from Manage reaction databases.

### Added. Engines & scoring
- **Rule-based retrosynthesis (RetroRules)**. `retrorules.py`: RDKit rule application,
  bounded retro-expansion grounding in the host's InChIKeys (with a wall-clock budget and
  an element pre-filter that cuts the ruleset ~12× for heteroatom targets), a bundled
  curated ruleset, and a full-dataset asset manager. `install_full_ruleset()` downloads,
  MD5-verifies and extracts the **full RetroRules RR02 dataset** (350k rules, 8
  diameters) and caches a compact per-diameter slice for fast reloads. The engine able
  to propose chemistry absent from every database; GUI wiring is the next step.

### Fixed. Escher (refined)
- Left-click now **selects** a reaction (for editing its flux in the FBA navigator) and
  never opens a popup; the details popup comes only from **right-click ▸ Show reaction
  details**. A left-click that moved (a pan/drag) is not treated as a selection.

## 0.2.1 - 2026-07-14

A correctness release for the **Pathway Design** engine. Five bugs are fixed that could
each make a designed pathway wrong, unusable, or simply missing. Also included are the
interface fixes for the Escher Visualizer and the main window.

### Fixed. Pathway Design engine
- **Cofactors written generically are now understood.** Databases such as ModelSEED
  write oxidation steps against an abstract `Reduced flavin`/`Flavin` pair rather than
  a concrete carrier. Cofactor detection only recognised BiGG-style ids, so "Reduced
  flavin" was mistaken for a compound that had to be synthesised from scratch and the
  reaction could never be used. This blocked every flavin-monooxygenase and P450 step. Including the oxidation of valencene to **nootkatol** and **nootkatone**, which are
  now found. When such a pathway is applied, the placeholder is mapped onto your model's
  real carrier (e.g. FADH2/FAD) and the substitution is reported, instead of entering the
  model as a dead-end compound that silently forces zero production.
- **Compounds can no longer appear from nothing.** The search is compartment-agnostic,
  so a transport reaction (`X_c <=> X_e`) collapsed to `X <=> X` and left an empty
  substrate requirement. Making *any* compound with a transport available in one step.
  Routes could therefore rely on an external supply of a compound the host cannot make
  (e.g. hexanoyl-CoA in Synechocystis). Choosing starting metabolites also silently
  enabled the database's own exchange reactions, the opposite of the intended narrowing;
  the starting list now restricts the search as documented.
- **The search is now deterministic.** The same query could previously return a route, a
  different route, or no route at all on successive runs.
- **Routes that cannot carry flux are explained.** The search checks that every
  substrate is reachable but not that by-products can be disposed of. Such dead ends are
  now detected and named in the result, and the engine retries around the offending step.
- **A precursor must be one the host can actually produce.** A metabolite present in the
  model may still be blocked (unable to carry any flux), and routes were bottoming out on
  such compounds. Precursors are now checked against the host's real production capacity.

### Fixed. Interface
- The window no longer grows beyond the screen or leaves its maximized state when an
  analysis runs or a strategy is saved.
- Escher: "Move nodes" now moves nodes; added a show/hide legend control; the map uses
  the full window width.
- Strategy Visualizer: the whole-model difference map is only computed when you click
  **Draw Network**; control labels are no longer cropped.

## 0.2.0 - 2026-07-14

The second release. Headline additions: an interactive **Escher Visualizer**, an
**Omics** workflow (dataset preparation + eFlux/GIMME), consortia objective controls,
a much richer Pathway Design flow, and a book-length illustrated manual.

### Added
- **Escher Visualizer** tab. Interactive, publication-quality metabolic maps powered
  by the bundled Escher engine (QtWebEngine). Auto-generates a clean, pathway-aware map
  for *any* model. Three modes: single-strategy flux, difference (A − B) with an
  interpretive up/down-regulation legend, and a live **FBA navigator** (click a reaction,
  edit its bounds or knock it out, and the model re-solves instantly). A clean custom
  toolbar (Fit / **Move nodes** / Export SVG / Export PNG), hover highlighting and
  click-for-details on reactions and metabolites.
- **Omics** tab and **Prepare omics dataset** tool. Parse transcriptomics / proteomics /
  metabolomics tables (CSV/TSV/Excel) and map their identifiers onto the model's genes,
  with a coverage summary. Feeds eFlux and GIMME for condition-specific predictions.
- **Consortia objective** settings. Per-member minimum-growth and dominance (abundance)
  weighting for community models, so a fast grower can't starve its partner.
- **Multi-select Focus** (Whole model / Category / Subsystem) on both the Escher
  Visualizer and the Network Map, for clean, focused views.
- General **product-export** options when applying a heterologous pathway (secreted /
  volatile gas / intracellular accumulation / none), and editable ids & names for the
  new reactions/metabolites right in the Add-Pathway dialog.
- Book-length **user manual** (docs/GSM_ToolBox_Manual.pdf, ~60 pages) with two detailed,
  fully worked example chapters: a CO₂-to-styrene *Synechocystis* + *P. Putida*
  consortium, and an omics-integration case study.

### Changed
- **Tools** is now a top-level menu (File · Edit · View · Settings · Tools · Help) instead
  of a tab. The visualization tabs were renamed **Strategy Visualizer** and **Escher
  Visualizer**.
- The Information panel is hidden by default (a Model Summary pops up on load; reopen it
  from *View ▸ Information*), keeping the window within the screen.
- Analysis method buttons wrap text to width; the Scope/Objective controls sit at the top.

### Fixed
- SBML `LPAREN/RPAREN` compartment encoding no longer leaks into the Explorer tables or
  reaction equations.
- Pathway Design "Draw Scheme" now linearises the route correctly (cofactors are no longer
  mistaken for pathway intermediates).
- Applying a pathway always guarantees a single, mass-balanced product exchange.
- The window no longer exceeds the screen or drops out of Maximized when running analyses
  or switching tabs.

## 0.1.0

Initial release: model loading/editing, the constraint-based analysis suite
(FBA/pFBA/FVA, envelopes, phase planes, deletions, strain design, FSEOF, gap-filling,
omics eFlux/GIMME, EFM, community modelling), heterologous Pathway Design against
BiGG/MetaNetX/ModelSEED/KEGG, the Network Map and Strategy Explorer, and a Windows
installer.
