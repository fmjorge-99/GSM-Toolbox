# GSM ToolBox

A desktop application for genome-scale metabolic modelling, written for biologists rather
than for programmers. It wraps [COBRApy](https://opencobra.github.io/cobrapy/) in an
interface where every analysis is a button, every result is a table you can read, and the
assumptions behind a number are stated next to it.

Runs on Windows, macOS and Linux.

## What it does

**Simulate.** FBA, pFBA, FVA, gene and reaction essentiality, phenotype phase planes,
production envelopes and robustness analysis.

**Design strains.** OptKnock, RobustKnock, minimal cut sets and FSEOF. Results are ranked
by carbon yield rather than raw flux, so routes to different products stay comparable.

**Design pathways.** Search BiGG, MetaNetX, ModelSEED, KEGG and Rhea for a route from your
host to a target compound. When no database has one, the toolbox can generate candidate
chemistry from reaction rules. Every proposed route is checked for mass and charge
balance, for thermodynamic feasibility, and for steps that quietly swap a compound for a
branched-chain isomer of it. That last kind of step balances perfectly and is still
impossible.

**Find enzymes.** Candidate enzymes by EC number, or by reaction similarity for
rule-generated steps that carry no EC number at all.

**Follow dynamics.** Sweep a condition to find where the network rewires, or follow a
batch culture as it consumes its medium. Runs are stored side by side and can be overlaid
on one plot.

**Apply regulation.** A stoichiometric model will fix carbon in the dark and keep
photosynthesising after nitrogen has run out. Loadable rule sets constrain it with what
the organism actually does. Before you read any result, the interface tells you how much
of a rule set can affect your model, because a rule can report as firing while changing
nothing.

**Visualise.** An interactive network map with flux overlays, Escher maps and
publication-ready figures.

## Install

| Your platform | Do this |
|---|---|
| **Windows** | Download the installer from [Releases](../../releases) and double-click it. Nothing else needed. |
| **macOS (Apple Silicon)** | Download the `arm64` `.dmg` from [Releases](../../releases), or install from source — see below. |
| **macOS (Intel)** | Install from source. (Described below) |
| **Linux** | Download the `.tar.gz` from [Releases](../../releases), or install from source — see below. |

### macOS

Installing from source takes about five minutes, works on Apple Silicon and Intel
alike, and is the recommended route on an Intel Mac — Apple has dropped that
architecture, so a prebuilt `x86_64` `.dmg` may not be in every release.

**1. Install Python 3.10, 3.11 or 3.12.** macOS ships 3.9, which is too old. Download the
macOS universal2 installer from [python.org](https://www.python.org/downloads/macos/) and
double-click it, or run `brew install "python@3.12"`.

**2. In Terminal:**

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
bash scripts/install.sh
```

That creates a double-clickable **`GSM ToolBox.app`** in the folder — drag it to
Applications. The script checks that the libraries really import and that Qt can start
before telling you it worked, so a broken install reports itself instead of failing at
first launch. Apple Silicon and Intel are both supported.

> Use `bash scripts/install.sh`, not `./scripts/install.sh`. The scripts are stored
> without the executable bit, so the `./` form gives `permission denied`.

**Want a self-contained app instead?** `bash scripts/build_bundle.sh` freezes Python, Qt
and the solvers into `dist/GSM ToolBox.app` plus a `.dmg` you can share with people who
have no Python. It needs Python once, to build with, and takes about fifteen minutes.

### Linux

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
bash scripts/install.sh
```

The script prints the exact Qt system packages your distribution needs if any are
missing, then adds an entry to your applications menu.

### Any platform, manually

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_gsm_toolbox.py
```

Full instructions are in [INSTALL.md](INSTALL.md): the system libraries Qt needs on
Linux, the Gatekeeper step for a downloaded macOS build, and how to build the Windows
installer, the macOS `.dmg` and the Linux tarball yourself.

## Documentation

| Document | What it covers |
|---|---|
| [`docs/GSM_ToolBox_Manual.pdf`](docs/GSM_ToolBox_Manual.pdf) | The user manual. Concepts, every analysis, worked examples |
| [`docs/pathway_manual/`](docs/pathway_manual) | Heterologous pathway design in depth |
| [`CHANGELOG.md`](CHANGELOG.md) | Release history |

## Models and data

No metabolic models are distributed with the toolbox. It is software, and the models
belong to their publishers.

Any SBML or JSON model COBRApy can read will open, including anything from
[BiGG Models](http://bigg.ucsd.edu/), the [BioModels](https://www.ebi.ac.uk/biomodels/)
repository, or your own reconstruction. A small example model is bundled so the interface
has something to open on a first run.

A compact universal reaction database ships with the application so pathway search works
without a network connection. Larger databases are downloaded from their providers when
you ask for them, and are not redistributed here.

## A note on what the numbers mean

Constraint-based models predict what a network can do at steady state, not what a cell
will do. The application is built around that distinction. Thresholds carry a stated
confidence, an inferred result says it is inferred, a measurement that is missing is
reported as missing rather than as zero, and a route whose chemistry looks wrong is
flagged even when it balances. Predictions are hypotheses to test at the bench.

## Licence

GSM ToolBox was designed and developed by **Jorge Fernández Méndez**. It is released
under GPL-3.0-or-later. See [`LICENSE`](LICENSE) for the terms and
[`AUTHORS.md`](AUTHORS.md) for authorship.

If you redistribute this software or build on it, the licence requires you to keep the
author notice and to make your source available under the same terms. Naming the original
author is part of that.

The pre-built binaries bundle solver and platform runtime libraries whose terms are not
GPL compatible. A linking exception under section 7 of the GPL covers that case and is
explained in [`LICENSE-EXCEPTION.md`](LICENSE-EXCEPTION.md). Every third-party licence is
recorded in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md). Installing from source
is unaffected.

## Citing

If this contributes to published work, please cite GSM ToolBox and its author, along with
the methods and data you used. That means COBRApy, the source publication of whichever
model you started from, and the providers of any database you fetched.
`THIRD_PARTY_LICENSES.md` lists them.
