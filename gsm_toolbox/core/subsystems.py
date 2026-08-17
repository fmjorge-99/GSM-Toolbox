"""Detect which central pathway each reaction belongs to, for models that say nothing.

Subsystem annotations are optional in SBML and frequently absent — reconstructions built
by automated pipelines, models converted between formats, and anything assembled from a
universal database often arrive with every ``subsystem`` field blank. That costs more than
tidiness: subsystem is what the Explorer groups by, what Escher focuses on, what FSEOF
ranks by, and what a regulatory rule targets when the user wants to affect "the CCM"
rather than eleven reaction ids.

**How detection works, and why it is layered.** There is no single reliable signal, so
five are used and each assignment records which one fired:

| Evidence | Strength | Why |
|---|---|---|
| ``structure`` | highest | A boundary reaction *is* an exchange; a reaction moving one compound between compartments without changing it *is* a transport. These are facts about the stoichiometry, not guesses. |
| ``id`` | high | BiGG ids for central metabolism are near-universal (``PGI``, ``GAPD``, ``RBPC``). An exact match is strong evidence. |
| ``ec`` | high | EC numbers are namespace-independent, so they carry across BiGG, MetaNetX, SEED and KEGG models. |
| ``metabolite`` | moderate | Presence of a pathway-specific intermediate (``6pgc`` for the OPP, ``skm`` for shikimate). Intermediates are specific; currency metabolites are excluded. |
| ``name`` | weakest | Keyword in the reaction name. Last resort, and flagged as such. |

Evidence is reported per reaction rather than hidden, because the user is expected to
review the result before it is written. A "glycolysis" assignment inferred from a name
containing "kinase" deserves a different level of scrutiny than one matched on EC 5.3.1.9,
and collapsing both into an unqualified answer would remove the reader's ability to tell.

**Nothing is written until the caller says so.** :func:`detect` is pure; :func:`apply`
performs the mutation. Existing subsystem annotations are preserved by default — a model
that already carries curated subsystems should not have them overwritten by inference.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import cobra

# --------------------------------------------------------------------------------------
# Evidence
# --------------------------------------------------------------------------------------
STRUCTURE = "structure"
BY_ID = "id"
BY_EC = "ec"
BY_METABOLITE = "metabolite"
BY_NAME = "name"

#: Strongest first. Used to resolve a reaction matching several pathways.
EVIDENCE_ORDER = [STRUCTURE, BY_ID, BY_EC, BY_METABOLITE, BY_NAME]

EVIDENCE_LABEL = {
    STRUCTURE: "stoichiometry",
    BY_ID: "reaction id",
    BY_EC: "EC number",
    BY_METABOLITE: "pathway intermediate",
    BY_NAME: "reaction name",
}

#: Metabolites that appear across the whole network and therefore identify nothing.
#: Without this list a metabolite signature would drag half the model into whichever
#: pathway happened to claim ATP.
CURRENCY = {
    "h", "h2o", "atp", "adp", "amp", "pi", "ppi", "nad", "nadh", "nadp", "nadph",
    "co2", "o2", "nh4", "coa", "accoa", "q8", "q8h2", "fad", "fadh2", "gtp", "gdp",
    "so4", "hco3", "na1", "k", "mg2", "fe2", "fe3", "cl", "h2o2", "no3", "no2",
    # S-adenosyl-methionine and its product are group-transfer cofactors, not markers
    # of one-carbon metabolism: measured against iJN678 they appear in porphyrin,
    # quinone and lipid methylations throughout the model.
    "amet", "ahcys",
}

#: The 20 proteinogenic amino acids, BiGG stems.
AMINO_ACIDS = {
    "ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L", "glu__L", "gly",
    "his__L", "ile__L", "leu__L", "lys__L", "met__L", "phe__L", "pro__L", "ser__L",
    "thr__L", "trp__L", "tyr__L", "val__L",
}

#: Nucleotide-pathway intermediates. Deliberately excludes ATP/GTP and the other
#: currency triphosphates, which carry no pathway information.
NUCLEOTIDE_INTERMEDIATES = {
    "prpp", "imp", "xmp", "gmp", "ump", "omp", "cmp", "dump", "dtmp", "dgmp", "damp",
    "dcmp", "air", "aicar", "faicar", "gar", "fgam", "fgam", "cair", "saicar", "pram",
    "dhor__S", "orot", "orot5p", "adn", "gsn", "ins", "ura", "thym", "csn", "din",
}


# --------------------------------------------------------------------------------------
# The catalogue
# --------------------------------------------------------------------------------------
@dataclass
class Pathway:
    """One central pathway and the signals that identify its reactions."""

    name: str
    description: str
    ids: Set[str] = field(default_factory=set)
    ecs: Set[str] = field(default_factory=set)
    metabolites: Set[str] = field(default_factory=set)
    keywords: Tuple[str, ...] = ()
    #: How many of this pathway's own intermediates a reaction must touch before the
    #: metabolite signature counts. Two is the useful default: one shared intermediate is
    #: usually a branch point into a different pathway.
    metabolite_hits: int = 2


#: Pathways expected in almost any organism, in the order they are presented.
#:
#: EC numbers ending in a dot are prefixes (``1.7.1.`` matches every nitrate/nitrite
#: reductase). Reaction ids are BiGG; a model in another namespace falls back to EC,
#: intermediates and finally names.
CATALOGUE: List[Pathway] = [
    Pathway(
        name="Glycolysis / Gluconeogenesis",
        description="Glucose ⇄ pyruvate, and the gluconeogenic reverse.",
        ids={"HEX1", "GLCpts", "PGI", "PFK", "FBP", "FBA", "TPI", "GAPD", "PGK",
             "PGM", "ENO", "PYK", "PPS", "PDH", "G3PD1", "FBA3", "PFK_3", "HEX7"},
        ecs={"2.7.1.1", "2.7.1.2", "5.3.1.9", "2.7.1.11", "3.1.3.11", "4.1.2.13",
             "5.3.1.1", "1.2.1.12", "2.7.2.3", "5.4.2.11", "5.4.2.12", "4.2.1.11",
             "2.7.1.40", "2.7.9.2", "5.4.2.1"},
        metabolites={"g6p", "f6p", "fdp", "dhap", "g3p", "13dpg", "3pg", "2pg", "pep",
                     "pyr", "glc__D"},
        keywords=("glycolys", "gluconeogen"),
    ),
    Pathway(
        name="Pentose phosphate pathway",
        description="Oxidative and non-oxidative branches — NADPH and R5P supply.",
        ids={"G6PDH2r", "PGL", "GND", "RPE", "RPI", "TKT1", "TKT2", "TALA", "ZWF",
             "G6PBDH", "PGLc"},
        ecs={"1.1.1.49", "1.1.1.363", "3.1.1.31", "1.1.1.44", "1.1.1.343", "5.1.3.1",
             "5.3.1.6", "2.2.1.1", "2.2.1.2"},
        metabolites={"6pgl", "6pgc", "ru5p__D", "xu5p__D", "r5p", "s7p", "e4p"},
        keywords=("pentose phosphate", "oxidative pp"),
    ),
    Pathway(
        name="Entner–Doudoroff pathway",
        description="The KDPG route from 6-phosphogluconate to pyruvate + G3P.",
        ids={"EDD", "EDA"},
        ecs={"4.2.1.12", "4.1.2.14"},
        metabolites={"2ddg6p", "6pgc"},
        keywords=("entner", "doudoroff"),
    ),
    Pathway(
        name="TCA cycle",
        description="Citrate cycle — oxidation of acetyl-CoA and precursor supply.",
        ids={"CS", "ACONT", "ACONTa", "ACONTb", "ICDHyr", "ICDHx", "AKGDH", "AKGDb",
             "SUCOAS", "SUCDi", "SUCD1", "FUM", "MDH", "OOR2r", "FRD7"},
        ecs={"2.3.3.1", "2.3.3.8", "4.2.1.3", "1.1.1.42", "1.1.1.41", "1.2.4.2",
             "2.3.1.61", "6.2.1.5", "1.3.5.1", "4.2.1.2", "1.1.1.37"},
        metabolites={"cit", "icit", "akg", "succoa", "succ", "fum", "mal__L", "oaa"},
        keywords=("tca", "citrate cycle", "citric acid cycle", "krebs"),
    ),
    Pathway(
        name="Glyoxylate shunt",
        description="Isocitrate lyase and malate synthase — the anaplerotic bypass.",
        ids={"ICL", "MALS"},
        ecs={"4.1.3.1", "2.3.3.9"},
        metabolites={"glx", "icit"},
        keywords=("glyoxylate",),
    ),
    Pathway(
        name="Anaplerosis and pyruvate metabolism",
        description="Carboxylation and decarboxylation reactions replenishing the TCA.",
        ids={"PPC", "PPCK", "PC", "ME1", "ME2", "MDH2", "ACKr", "PTAr", "LDH_D",
             "ALCD2x", "ACALD", "POR5", "PFL"},
        ecs={"4.1.1.31", "4.1.1.49", "6.4.1.1", "1.1.1.38", "1.1.1.40", "2.7.2.1",
             "2.3.1.8", "1.1.1.28", "1.1.1.1", "1.2.1.10", "2.3.1.54"},
        metabolites={"pyr", "oaa", "actp", "ac", "lac__D", "etoh"},
        keywords=("anaplero", "pyruvate metabolism"),
    ),
    Pathway(
        name="Carbon fixation (Calvin–Benson cycle)",
        description="RuBisCO, phosphoribulokinase and the regenerative phase.",
        ids={"RBPC", "RBCh", "PRUK", "SBP", "SBTA", "FBA3", "GAPDH_nadp", "TKT1",
             "TKT2", "RPE", "RPI", "PGK"},
        ecs={"4.1.1.39", "2.7.1.19", "3.1.3.37", "1.2.1.13", "3.1.3.11"},
        metabolites={"rb15bp", "ru5p__D", "sbp", "s17bp", "3pg", "e4p"},
        keywords=("calvin", "carbon fixation", "rubisco", "cbb"),
    ),
    Pathway(
        name="Photosynthesis and light reactions",
        description="Photosystems, electron transport and photophosphorylation.",
        ids={"PSI", "PSII", "PSIIh", "CBFC", "CBFCumq", "FNOR", "CEF", "ATPSh",
             "NDH1", "NDH2", "PQH2t", "FQR", "CYO1b"},
        ecs={"1.10.3.9", "1.18.1.2", "7.1.1.-"},
        metabolites={"photon", "pq", "pqh2", "pc", "pcox", "fdxo", "fdxr"},
        keywords=("photosystem", "photosynth", "light reaction", "plastoquinone",
                  "photophosphoryl"),
    ),
    Pathway(
        name="Oxidative phosphorylation",
        description="Respiratory chain and the ATP synthase.",
        ids={"NADH16", "NADH17", "NADH18", "CYTBD", "CYTBO3_4pp", "ATPS4r", "ATPS4rpp",
             "CYOO", "CYOR", "THD2", "NADTRHD", "ATPM"},
        ecs={"7.1.1.2", "7.1.1.7", "7.1.1.9", "7.1.2.2", "3.6.3.14", "1.6.1.2"},
        metabolites={"q8", "q8h2", "mql8", "mqn8"},
        keywords=("oxidative phosphoryl", "atp synthase", "respirat",
                  "electron transport"),
    ),
    Pathway(
        name="Nitrogen assimilation",
        description="Nitrate/nitrite/ammonium uptake into glutamate and glutamine.",
        ids={"GLNS", "GLUDy", "GLUDx", "GLUSy", "GLUSx", "NO3R1", "NO3R2", "NTRIR2x",
             "NTRIR4pp", "NIT", "UREA", "URCN", "UREASE", "ASNS1", "ASNS2", "N2Ored"},
        ecs={"6.3.1.2", "1.4.1.13", "1.4.1.14", "1.4.1.2", "1.4.1.3", "1.4.1.4",
             "1.7.1.", "1.7.7.1", "1.7.7.2", "3.5.1.5", "1.18.6.1", "6.3.5.4"},
        metabolites={"no3", "no2", "nh4", "gln__L", "glu__L", "urea", "cbp"},
        keywords=("nitrogen assimil", "nitrate", "nitrite reduct", "glutamine synth",
                  "nitrogen metabolism", "urease"),
    ),
    Pathway(
        name="Sulfur assimilation",
        description="Sulfate reduction to sulfide and cysteine formation.",
        ids={"SULabc", "SADT2", "ADSK", "PAPSR", "BPNT", "SULR", "CYSS", "SERAT"},
        ecs={"2.7.7.4", "2.7.1.25", "1.8.4.8", "1.8.1.2", "2.5.1.47", "2.3.1.30"},
        metabolites={"so4", "aps", "paps", "so3", "h2s", "acser"},
        keywords=("sulfate", "sulfur", "sulphate"),
    ),
    Pathway(
        name="Shikimate pathway",
        description="E4P + PEP to chorismate — the aromatic precursor route.",
        ids={"DDPA", "DHQS", "DHQTi", "DHQD", "SHK3Dr", "SHKK", "PSCVT", "CHORS",
             "CHORM", "PPND", "PPNDH"},
        ecs={"2.5.1.54", "4.2.3.4", "4.2.1.10", "1.1.1.25", "2.7.1.71", "2.5.1.19",
             "4.2.3.5"},
        metabolites={"2dda7p", "3dhq", "3dhsk", "skm", "skm5p", "3psme", "chor"},
        keywords=("shikimate", "chorismate"),
    ),
    Pathway(
        name="Amino acid metabolism",
        description="Biosynthesis and degradation of the proteinogenic amino acids.",
        ids={"ASPTA", "ALATA_L", "VALTA", "ILETA", "LEUTA", "PHETA1", "TYRTA",
             "GLUN", "ASPK", "ASAD", "HSDy", "HSK", "THRS", "DAPDC", "DAPE",
             "ACLS", "KARA1", "DHAD1", "IPPS", "IPPMIa", "OMCDC", "ANS", "ANPRT",
             "PRAIi", "IGPS", "TRPS1", "TRPS2", "TRPS3", "HISTD", "HISTP", "IG3PS",
             "ACGS", "ACGK", "AGPR", "ACOTA", "ORNTA", "OCBT", "ARGSS", "ARGSL",
             "SDPTA", "SHSL1", "METS", "CYSTL", "PROD2", "P5CR", "G5SD", "GLU5K"},
        ecs={"2.6.1.1", "2.6.1.2", "2.6.1.42", "2.6.1.6", "2.7.2.4", "1.2.1.11",
             "1.1.1.3", "2.7.1.39", "4.2.3.1", "4.1.1.20", "2.2.1.6", "1.1.1.86",
             "4.2.1.9", "2.3.3.13", "4.2.1.33", "1.1.1.85", "4.1.3.27", "2.4.2.18",
             "5.3.1.24", "4.1.1.48", "4.2.1.20", "2.6.1.11", "2.1.3.3", "6.3.4.5",
             "4.3.2.1", "2.5.1.48", "4.4.1.8", "2.1.1.14", "1.5.1.2", "1.2.1.41"},
        metabolites=set(AMINO_ACIDS),
        keywords=("amino acid", "aminotransferase", "transaminase"),
        metabolite_hits=2,
    ),
    Pathway(
        name="Nucleotide biosynthesis (DNA/RNA)",
        description="Purine and pyrimidine synthesis, salvage and deoxyribonucleotides.",
        ids={"PRPPS", "GLUPRT", "PRAGSr", "GARFT", "PRFGS", "PRAIS", "AIRC2",
             "AIRC3", "PRASCS", "ADSL2r", "AICART", "IMPC", "ADSS", "ADSL1r",
             "IMPD", "GMPS2", "ASPCT", "DHORTS", "DHORD2", "ORPT", "OMPDC", "UMPK",
             "NDPK1", "NDPK2", "NDPK3", "NDPK4", "CTPS2", "RNDR1", "RNDR2", "RNDR3",
             "RNDR4", "TMDS", "DTMPK", "URIDK2r", "ADK1", "GK1", "CYTK1", "NTD",
             "PUNP1", "HXPRT", "GUAPRT", "ADPT", "UPPRT"},
        ecs={"2.7.6.1", "2.4.2.14", "6.3.4.13", "2.1.2.2", "6.3.5.3", "6.3.3.1",
             "4.1.1.21", "6.3.2.6", "4.3.2.2", "2.1.2.3", "3.5.4.10", "6.3.4.4",
             "1.1.1.205", "6.3.5.2", "2.1.3.2", "3.5.2.3", "1.3.5.2", "2.4.2.10",
             "4.1.1.23", "2.7.4.22", "2.7.4.6", "6.3.4.2", "1.17.4.1", "2.1.1.45",
             "2.7.4.9", "2.7.4.3", "2.7.4.8", "2.7.4.14", "2.4.2.8", "2.4.2.7",
             "2.4.2.22", "2.4.2.9"},
        metabolites=set(NUCLEOTIDE_INTERMEDIATES),
        keywords=("purine", "pyrimidine", "nucleotide", "ribonucleotide",
                  "deoxyribonucle", "thymidylate"),
    ),
    Pathway(
        name="Carbon storage (glycogen / starch / PHB)",
        description="Polymeric carbon reserves and their mobilisation.",
        ids={"GLGC", "GLCS1", "GLBRAN2", "GLDBRAN2", "GLPASE1", "GLPASE2", "GLCP",
             "AMALT1", "PGMT", "PHAS", "ACACT1r", "PHAR", "AGPAT"},
        # No "2.3.1.-" wildcard: measured against iJN678 it pulled in 37 fatty-acid
        # synthase reactions (EC 2.3.1.41/85/86), which are a different pathway entirely.
        ecs={"2.7.7.27", "2.4.1.21", "2.4.1.18", "2.4.1.1", "5.4.2.2", "2.3.1.9",
             "1.1.1.36"},
        metabolites={"glycogen", "glgc", "starch", "g1p", "malt", "phb", "3hbcoa"},
        keywords=("glycogen", "starch", "polyhydroxy", "carbon storage", "granule"),
    ),
    Pathway(
        name="Fatty acid and lipid metabolism",
        description="Fatty-acid synthesis and elongation, and membrane lipids.",
        ids={"ACCOAC", "MCOATA", "KAS14", "KAS15", "3OAR40", "3HAD40", "EAR40x",
             "FACOAL140", "AGPAT160", "PSSA120", "PSD120", "CDAPPA160", "G3PAT120"},
        ecs={"6.4.1.2", "2.3.1.39", "2.3.1.41", "1.1.1.100", "4.2.1.59", "1.3.1.9",
             "6.2.1.3", "2.3.1.15", "2.3.1.51", "2.7.8.5", "2.7.8.8", "4.1.1.65"},
        metabolites={"malacp", "acacp", "accoa", "3hdecacp", "ddca", "ttdca", "hdca"},
        keywords=("fatty acid", "lipid", "phospholipid", "acyl carrier"),
    ),
    Pathway(
        name="Cell wall and membrane biogenesis",
        description="Peptidoglycan, lipopolysaccharide and murein assembly.",
        ids={"UAGDP", "UAGCVT", "UAPGR", "UAMAS", "UAMAGS", "UGMDDS", "UAAGDS",
             "PAPPT3", "MPTG", "MURI", "GLUR", "ALAALAr", "ALARi"},
        ecs={"2.7.7.23", "2.5.1.7", "1.3.1.98", "6.3.2.8", "6.3.2.9", "6.3.2.10",
             "6.3.2.13", "2.7.8.13", "2.4.1.129", "5.1.1.3", "5.1.1.1", "6.3.2.4"},
        metabolites={"uacgam", "uamr", "uama", "uamag", "murein", "lps", "kdo2lipid4"},
        keywords=("peptidoglycan", "murein", "cell wall", "lipopolysaccharide",
                  "cell envelope"),
    ),
    Pathway(
        name="Cofactor and vitamin biosynthesis",
        description="NAD, CoA, folate, thiamine, riboflavin, haem and quinones.",
        ids={"NNDPR", "NNAT", "NADS1", "NADK", "PANTS", "PPNCL2", "PPCDC", "DPCOAK",
             "DHFR", "DHFS", "DHPS2", "GTPCI", "THZPSN", "TMPPP", "RBFSa", "RBFSb",
             "RIBFLVt", "HEMEOS", "CPPPGO", "PPBNGS", "UPP3S", "SHCHF", "ADCS"},
        ecs={"2.4.2.19", "2.7.7.18", "6.3.1.5", "2.7.1.23", "6.3.2.1", "4.1.1.36",
             "2.7.1.24", "1.5.1.3", "6.3.2.12", "2.5.1.15", "3.5.4.16", "2.5.1.3",
             "2.7.6.2", "2.5.1.9", "3.5.4.25", "4.2.3.5", "2.5.1.19", "4.2.1.75",
             "4.1.1.37", "2.5.1.61", "1.3.3.3", "4.99.1.1"},
        metabolites={"nad", "nmn", "dhf", "thf", "10fthf", "5mthf", "pan4p", "dpcoa",
                     "ribflv", "fmn", "thmpp", "ppbng", "uppg3", "pheme", "sheme"},
        keywords=("cofactor", "vitamin", "folate", "thiamin", "riboflavin", "biotin",
                  "haem", "heme", "porphyrin", "quinone biosynth", "pantothenate",
                  "coenzyme a biosynth"),
    ),
    Pathway(
        name="One-carbon and methyl-group metabolism",
        description="Folate-mediated C1 transfer and the methionine cycle.",
        ids={"GHMT2r", "MTHFC", "MTHFD", "MTHFR2", "FTHFD", "METAT", "AHCi", "SAHH",
             "FOMETRi"},
        ecs={"2.1.2.1", "3.5.4.9", "1.5.1.5", "1.5.1.20", "3.5.1.10", "2.5.1.6",
             "3.3.1.1", "2.1.1.13"},
        metabolites={"methf", "mlthf", "5mthf", "10fthf", "amet", "ahcys", "hcys__L"},
        keywords=("one-carbon", "one carbon", "methionine cycle"),
    ),
    Pathway(
        name="Transport",
        description="Movement of a compound between compartments.",
        keywords=("transport", "permease", "symport", "antiport", "abc", "diffusion",
                  "uniport", "efflux", "uptake"),
    ),
    Pathway(
        name="Exchange / demand / sink",
        description="Boundary reactions defining what crosses the system boundary.",
        keywords=("exchange", "demand", "sink"),
    ),
    Pathway(
        name="Biomass and maintenance",
        description="Growth objective and maintenance requirements.",
        keywords=("biomass", "growth", "maintenance"),
    ),
]

_PATHWAY_BY_NAME = {p.name: p for p in CATALOGUE}


# --------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------
@dataclass
class Assignment:
    reaction_id: str
    reaction_name: str
    pathway: str
    evidence: str                # one of the evidence constants
    detail: str = ""             # what actually matched, e.g. "EC 5.3.1.9"

    @property
    def strength(self) -> int:
        try:
            return EVIDENCE_ORDER.index(self.evidence)
        except ValueError:
            return len(EVIDENCE_ORDER)

    def evidence_text(self) -> str:
        label = EVIDENCE_LABEL.get(self.evidence, self.evidence)
        return f"{label}: {self.detail}" if self.detail else label


@dataclass
class DetectionReport:
    assignments: Dict[str, Assignment] = field(default_factory=dict)
    unassigned: List[str] = field(default_factory=list)
    #: Reactions that already carried a subsystem, left alone unless overwrite was asked.
    preserved: Dict[str, str] = field(default_factory=dict)
    total_reactions: int = 0

    def by_pathway(self) -> Dict[str, List[Assignment]]:
        out: Dict[str, List[Assignment]] = {}
        for assignment in self.assignments.values():
            out.setdefault(assignment.pathway, []).append(assignment)
        for items in out.values():
            items.sort(key=lambda a: a.reaction_id)
        return dict(sorted(out.items(), key=lambda kv: -len(kv[1])))

    def coverage(self) -> float:
        if not self.total_reactions:
            return 0.0
        covered = len(self.assignments) + len(self.preserved)
        return covered / self.total_reactions

    def summary(self) -> str:
        pathways = self.by_pathway()
        text = (f"{len(self.assignments)} of {self.total_reactions} reactions assigned "
                f"to {len(pathways)} pathway(s).")
        if self.preserved:
            text += (f" {len(self.preserved)} already had a subsystem and were left "
                     f"unchanged.")
        if self.unassigned:
            text += f" {len(self.unassigned)} could not be placed."
        weak = sum(1 for a in self.assignments.values() if a.evidence == BY_NAME)
        if weak:
            text += (f" {weak} rest on a name keyword only — the weakest evidence, "
                     f"worth reviewing.")
        return text


# --------------------------------------------------------------------------------------
# Detection
# --------------------------------------------------------------------------------------
def _stem(metabolite_id: str) -> str:
    """Strip the compartment suffix: ``g6p_c`` → ``g6p``."""
    return re.sub(r"_[a-z][a-z0-9]?$", "", metabolite_id)


def _base_id(reaction_id: str) -> str:
    """Undo the decorations enzyme-constrained models add to reaction ids."""
    text = re.sub(r"_TG_(forward|reverse)$", "", reaction_id)
    text = re.sub(r"(_copy\d+|_\d+)$", "", text)
    return text


def _ec_matches(ecs: Iterable[str], wanted: Set[str]) -> Optional[str]:
    for ec in ecs:
        if ec in wanted:
            return ec
        for candidate in wanted:
            # A trailing dot marks a prefix: "1.7.1." matches "1.7.1.4".
            if candidate.endswith(".") and ec.startswith(candidate):
                return ec
            if candidate.endswith(".-") and ec.startswith(candidate[:-1]):
                return ec
    return None


def _is_transport(rxn: cobra.Reaction) -> bool:
    """A reaction moving a compound between compartments without changing it.

    Structural, so it holds regardless of naming convention: the same chemical stem
    appears on both sides, and more than one compartment is involved.
    """
    if rxn.boundary or len(rxn.metabolites) < 2:
        return False
    compartments = {m.compartment for m in rxn.metabolites}
    if len(compartments) < 2:
        return False
    consumed = {_stem(m.id) for m, c in rxn.metabolites.items() if c < 0}
    produced = {_stem(m.id) for m, c in rxn.metabolites.items() if c > 0}
    shared = (consumed & produced) - CURRENCY
    return bool(shared)


def _looks_like_biomass(rxn: cobra.Reaction) -> bool:
    text = f"{rxn.id} {rxn.name or ''}".upper()
    return "BIOMASS" in text or rxn.id.upper() in {"ATPM", "NGAM"}


def _structural_pathway(rxn: cobra.Reaction) -> Optional[Tuple[str, str]]:
    """Pathways decidable from stoichiometry alone. Highest confidence available."""
    if _looks_like_biomass(rxn):
        return "Biomass and maintenance", "objective or maintenance reaction"
    if rxn.boundary:
        return "Exchange / demand / sink", "boundary reaction"
    if _is_transport(rxn):
        return "Transport", "same compound in two compartments"
    return None


def _candidates(rxn: cobra.Reaction, ecs: Sequence[str]) -> List[Assignment]:
    """Every pathway this reaction could belong to, with the evidence for each."""
    found: List[Assignment] = []
    base = _base_id(rxn.id)
    name = (rxn.name or "").lower()
    stems = {_stem(m.id) for m in rxn.metabolites} - CURRENCY

    for pathway in CATALOGUE:
        if pathway.ids and base in pathway.ids:
            found.append(Assignment(rxn.id, rxn.name or "", pathway.name, BY_ID, base))
            continue
        ec = _ec_matches(ecs, pathway.ecs) if pathway.ecs else None
        if ec:
            found.append(Assignment(rxn.id, rxn.name or "", pathway.name, BY_EC,
                                    f"EC {ec}"))
            continue
        if pathway.metabolites:
            hits = stems & pathway.metabolites
            if len(hits) >= pathway.metabolite_hits:
                found.append(Assignment(rxn.id, rxn.name or "", pathway.name,
                                        BY_METABOLITE,
                                        ", ".join(sorted(hits)[:3])))
                continue
        if pathway.keywords and name:
            hit = next((k for k in pathway.keywords if k in name), None)
            if hit:
                found.append(Assignment(rxn.id, rxn.name or "", pathway.name, BY_NAME,
                                        f"'{hit}'"))
    return found


def detect(model: cobra.Model, *, overwrite: bool = False,
           minimum_evidence: str = BY_NAME) -> DetectionReport:
    """Assign each reaction to a central pathway, without touching the model.

    ``overwrite=False`` (the default) leaves reactions that already declare a subsystem
    alone: a curated annotation is better evidence than anything inferred here, and
    replacing it would quietly discard the model author's work.

    ``minimum_evidence`` drops assignments weaker than the given level, for a caller that
    wants only high-confidence results.
    """
    from .databases import reaction_ec_numbers

    report = DetectionReport(total_reactions=len(model.reactions))
    try:
        floor = EVIDENCE_ORDER.index(minimum_evidence)
    except ValueError:
        floor = len(EVIDENCE_ORDER)

    for rxn in model.reactions:
        existing = (rxn.subsystem or "").strip()
        if existing and not overwrite:
            report.preserved[rxn.id] = existing
            continue

        structural = _structural_pathway(rxn)
        if structural:
            pathway, detail = structural
            report.assignments[rxn.id] = Assignment(
                rxn.id, rxn.name or "", pathway, STRUCTURE, detail)
            continue

        try:
            ecs = reaction_ec_numbers(rxn)
        except Exception:  # noqa: BLE001 — annotation shapes vary; never fail detection
            ecs = []
        options = [a for a in _candidates(rxn, ecs) if a.strength <= floor]
        if not options:
            report.unassigned.append(rxn.id)
            continue
        # Strongest evidence wins; ties break on catalogue order, which puts the more
        # specific pathway first (Entner–Doudoroff before glycolysis, for instance).
        options.sort(key=lambda a: (a.strength,
                                    [p.name for p in CATALOGUE].index(a.pathway)))
        report.assignments[rxn.id] = options[0]

    return report


def alternatives(model: cobra.Model, reaction_id: str) -> List[Assignment]:
    """Every pathway a reaction could plausibly belong to, for the review dialog."""
    from .databases import reaction_ec_numbers

    if not model.reactions.has_id(reaction_id):
        return []
    rxn = model.reactions.get_by_id(reaction_id)
    structural = _structural_pathway(rxn)
    out = []
    if structural:
        out.append(Assignment(rxn.id, rxn.name or "", structural[0], STRUCTURE,
                              structural[1]))
    try:
        ecs = reaction_ec_numbers(rxn)
    except Exception:  # noqa: BLE001
        ecs = []
    out.extend(_candidates(rxn, ecs))
    out.sort(key=lambda a: a.strength)
    return out


# --------------------------------------------------------------------------------------
# Applying
# --------------------------------------------------------------------------------------
def apply(model: cobra.Model, assignments: Dict[str, str]) -> int:
    """Write ``{reaction_id: subsystem}`` onto the model. Returns the number changed.

    Separate from :func:`detect` so nothing is written until the user has reviewed it.
    An empty string clears a subsystem, which is how the dialog expresses "leave this
    reaction unassigned".
    """
    changed = 0
    for rid, subsystem in assignments.items():
        if not model.reactions.has_id(rid):
            continue
        rxn = model.reactions.get_by_id(rid)
        new = (subsystem or "").strip()
        if (rxn.subsystem or "") != new:
            rxn.subsystem = new
            changed += 1
    return changed


def existing_subsystems(model: cobra.Model) -> Dict[str, int]:
    """Subsystems already present, with reaction counts."""
    counts: Dict[str, int] = {}
    for rxn in model.reactions:
        name = (rxn.subsystem or "").strip()
        if name:
            counts[name] = counts.get(name, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))
