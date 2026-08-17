# Third-party components

GSM ToolBox is GPL-3.0-or-later. This file records the licences of the libraries it
depends on, and of the additional binaries that appear inside the pre-built Windows
installer.

Versions are those the application was built and tested against. `pip` may resolve
newer ones; the licences below are the upstream projects' own and are not restated by
this project.

## Runtime dependencies (installed by `pip`)

| Component | Licence | GPL-3.0 compatible |
|---|---|---|
| PySide6 / shiboken6 (Qt for Python) | LGPL-3.0 | Yes |
| COBRApy | LGPL-2.0-or-later **or** GPL-2.0-or-later | Yes (via "or later") |
| optlang | Apache-2.0 | Yes (one-way, into GPL-3.0) |
| swiglpk / **GLPK** | **GPL-3.0** | Yes, and requires GPL compatibility |
| PySCIPOpt (binding) | MIT | Yes |
| SCIP (solver, version 9 or newer) | Apache-2.0 | Yes |
| StrainDesign | Apache-2.0 | Yes |
| RDKit | BSD-3-Clause | Yes |
| NumPy, SciPy, pandas, Pillow, lxml, cloudpickle | BSD-3-Clause | Yes |
| matplotlib | matplotlib licence (BSD-style) | Yes |
| python-libsbml | LGPL-2.1 | Yes |
| networkx | BSD-3-Clause | Yes |
| requests, depinfo | Apache-2.0 | Yes |
| openpyxl | MIT | Yes |
| equilibrator-api (optional) | MIT | Yes |

GLPK is the reason this project must be GPL compatible at all. It is GPL-3.0, and
COBRApy uses it as the default linear solver, so any distributed binary that links it
inherits that obligation. GPL-3.0-or-later is not simply a preference here. It is what
the dependency stack already requires.

## Additional binaries inside the pre-built Windows installer

These arrive inside the **PySCIPOpt** wheel and are bundled by PyInstaller. They are
*not* present in this repository.

| Binary | Component | Licence | GPL-3.0 compatible |
|---|---|---|---|
| `libscip.dll` | SCIP | Apache-2.0 | Yes |
| `coinmumps-3.dll` | MUMPS | Permissive, CeCILL-C style | Yes |
| `ipopt-3.dll` | Ipopt | **Eclipse Public License** | **No** |
| `libiomp5md.dll` | Intel OpenMP runtime | **Proprietary redistributable** | **No** |
| `libifcoremd.dll`, `libmmd.dll`, `svml_dispmd.dll` | Intel compiler runtimes | **Proprietary redistributable** | **No** |
| `msvcp140.dll` | Microsoft C++ runtime | Proprietary (system library) | Exempt as a System Library |

The last row is covered by the GPL's own **System Libraries** exception. The rows above
it are not, which is why [`LICENSE-EXCEPTION.md`](LICENSE-EXCEPTION.md) grants an
additional permission under GPL section 7 for the pre-built binary.

None of this affects a source installation: `pip` fetches those libraries onto your own
machine, which is not redistribution.

## Data

**No metabolic models are distributed with this software.** Two data files are bundled,
both functional parts of the application rather than research content:

| Asset | Origin | Terms |
|---|---|---|
| `gsm_toolbox/resources/databases/offline_universal.json` | Reaction database derived from the BiGG universal model | Cite [BiGG Models](http://bigg.ucsd.edu/) if you publish results obtained through it |
| `gsm_toolbox/resources/models/e_coli_core.xml` | The *E. coli* core model, the standard COBRApy demonstration model | Freely redistributable; cite Orth *et al.* 2010, *EcoSal Plus* |

The second exists only so *File ▸ Open Example Model* has something to open on a first
run. It is not required for any analysis.

Reaction databases fetched at runtime are downloaded from their own providers at your
request and are not redistributed here. These are BiGG, MetaNetX, ModelSEED, KEGG, Rhea,
RetroRules and UniProt. Each carries its own terms, and KEGG in
particular restricts automated bulk access and commercial use. Check them before relying
on any of these beyond personal research.
