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

Ready-made bundles are on the [Releases](../../releases) page. They contain Python, Qt and
the solvers, so nothing else needs installing.

| Platform | Get it |
|---|---|
| Windows, installer | Run `GSM_ToolBox_Setup_<version>.exe`. Installs per user, no admin rights |
| Windows, portable | Unzip `GSM_ToolBox-<version>-windows-portable.zip` and run `GSM_ToolBox.exe`. Nothing is installed |
| macOS | Open the `.dmg` and drag the app to Applications |
| Linux | `tar -xzf GSM_ToolBox-*-linux-x86_64.tar.gz` then `./GSM_ToolBox/GSM_ToolBox` |

To install from source on Linux or macOS, one command sets up the environment, checks that
Qt really works and adds a desktop launcher:

```sh
./scripts/install.sh && ./scripts/run.sh
```

Full instructions are in [INSTALL.md](INSTALL.md), including the system libraries Qt needs
on Linux and the extra step macOS requires on first launch.

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

Ideated and implemented by Jorge Fernández Méndez. Released under GPL-3.0-or-later. See
[`LICENSE`](LICENSE) for the terms and [`AUTHORS.md`](AUTHORS.md) for authorship.

If you redistribute this software or build on it, the licence requires you to keep the
copyright notice and to make your source available under the same terms. Attribution to
the original author is part of that.

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
