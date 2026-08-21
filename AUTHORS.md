# Author

GSM ToolBox was designed, developed and is maintained by **Jorge Fernández Méndez**.

The idea, the interface, the analysis workflows and the scientific decisions behind them
are his work. If the toolbox is useful to you, please credit him by name.

## Attribution

The licence asks that the author be named wherever this software, or work derived from
it, is redistributed or described. In practice that means keeping the author notice in
the source files, naming GSM ToolBox and its author in any publication that used it, and
keeping this file with the code.

## A note on the licence notice

The licence text uses the standard wording that the GNU General Public License expects,
which includes a copyright line naming the author. That line is what makes the licence
work: a person can only grant others permission to use, modify and share software if they
hold the rights to it in the first place. In most countries those rights arise
automatically when the work is created, with nothing to register and no fee to pay.

Two situations are worth checking before publishing more widely. If the software was
written as part of paid employment, or under the terms of a studentship or grant, your
employer or institution may hold or share those rights, and their policy decides what the
notice should say. This file is not legal advice, and an institutional research office
can answer both questions quickly.

## Contributing

Bug reports, feature requests and pull requests are welcome through the repository issue
tracker.

Contributions are accepted under the same licence as the project, GPL-3.0-or-later. See
[`LICENSE`](LICENSE). By opening a pull request you agree that your contribution may be
distributed under those terms. Please add yourself to this file in the same change.

## Acknowledgements

The application is built on open scientific software and does little that those projects
had not already made possible. COBRApy provides the constraint-based analysis, optlang
with GLPK and SCIP the optimisation, RDKit the chemistry, Qt the interface, and Escher
the pathway maps. Their licences, and those of every other dependency, are recorded in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

Reaction and compound data are fetched at the user's request from BiGG Models, MetaNetX,
ModelSEED, KEGG, Rhea, RetroRules and UniProt. None of it is redistributed here. Please
cite the providers whose data you use.
