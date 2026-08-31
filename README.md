# GSM Toolbox

**An interactive desktop application for genome-scale metabolic modeling, built for biologists.

[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3+-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)
![Platforms](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)

---

## What it is

GSM Toolbox is an interactive graphical environment for **loading, editing, constraining and simulating genome-scale metabolic models**. Every analysis, including customized ones, is run through buttons, selections, setting panels or guided wizards. The software is a toolbox that aims to unify the most common analysis for genome scale metabolic models under a single interface. 

**Intended for biologists with little to none coding skills**

**No command line is required, and no file is ever hand-written.** Configuration, constraint files and exports are produced by the interface itself.

Results are presented upfront, are navigable, and include plots and interactive visualizations that communicate better than numbers. However, raw data from every analysis is also easily accesible, exportable and reproducible. The aim is to **answer biologically meaningful questions** via custom analysis that are straightforward to set up, and deliver results that can be interpreted and acted on.

---

## Who it is for

Metabolic engineers, microbiologists and biotechnologists who want to use constraint-based modeling in their work, who do not write code, and who should not have to in order to run methods and bioinformatic tools that have been developed over the years, but always hard to access due to its scattered access and high entry barrier. The toolbox aim to be equally useful to people who *do* code, as a fast way to explore a model behaviour before committing to an script-based custom analysis. The **Export → Python script** functions is designed to hand you a working COBRApy script for any analysis that have been set up in the interface. Easing further implementations over the already generated code. 

---

## What it does

**Simulation**, flux balance analysis, parsimonious FBA, flux variability analysis, flux sampling, phenotype phase planes, production envelopes, robustness analysis, MOMA and ROOM, batch and dynamic simulation.

**Strain design**, OptKnock, RobustKnock, minimal cut sets and FSEOF, with results ranked by carbon yield so routes to different products are comparable.

**Pathway discovery**, searches BiGG, MetaNetX, ModelSEED, KEGG and Rhea for a route to a target compound, checking every candidate for mass and charge balance, thermodynamic feasibility, and spurious isomer swaps. Finds candidate enzymes by EC number or reaction similarity.

**Experimental data**, import measurements with their units and uncertainty, convert them to model constraints, and compare predictions against observations. Includes sealed-vessel gas handling and photon-flux conversion for phototrophs.

**Model quality and comparison**, MEMOTE quality reports, FROG reproducibility archives, energy-generating-cycle detection, curation tools, and the ability to run one analysis across several models at once and see where they agree and where they do not.

**Results and reproducibility**, every run is stored with full provenance, sessions save and reload everything, results can be re-run and verified, and anything can be exported as data, figures, a report, or a standalone Python script.

---


## AI Use Consideration

**GSM Toolbox code has been written using agentic AI**, specifically Claude Opus 5 and Sonnet 5, under human direction. The entire Toolbox interface design, workflows and and scientific judgment embodied in the toolbox has been carefully considere by a human. The AI agent acted as a translator from human language into code. Debugging, and final implementation has been always performed with an strict human-in-the-loop strategy for final curation.

During the last years the use of AI for solving scientifically relevant questions has exploded, and it is going to keep increasing. However te use of AI also carries serious considerations in terms of environmental impact, trustworthiness, reproducibility and intellectual property. All these aspects carry a significant ethical side, and the development of this toolbox has been done with all of them on mind.

The main inspiration for developing this toolbox is to facilitate access to Genome Scale Metabolic model analysis to users with little to none coding notions but deep knowledge about biological systems. Genome-scale metabolic analysis has well-established, deterministic methods that can be automated and audited. However, frequently code-naive users tend to directly use a language model to perform such an analysis. This not only significantly increases the use of AI but also produces results that cannot be easily reproduced or checked. GSM Toolbox exists so that a biologist without coding experience can run these analyses **properly**, through software whose behavior is fixed, inspectable and reproducible, rather than delegating them to a model that may improvise and consume significantly more resources. The software is tested against reference implementations: every wrapped analysis is asserted to return the same numbers as COBRApy, under more than one solver, in continuous integration. 

Using AI to build the tool, and then not needing AI to use it, is the intended trade.


---


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
