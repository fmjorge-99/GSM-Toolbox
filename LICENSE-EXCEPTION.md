# Additional permission under GNU GPL version 3 section 7

GSM ToolBox was designed and developed by Jorge Fernández Méndez, and is licensed under
the GNU General Public License, version 3 or, at your option, any later version. See
[`LICENSE`](LICENSE).

The following additional permission is granted under section 7 of that licence. It
applies to GSM ToolBox itself; it does not, and cannot, change the terms of any
third-party component.

> **Linking exception.** As a special exception, Jorge Fernández Méndez, as author and
> rights holder of GSM ToolBox, gives you permission to combine GSM ToolBox with the numerical solver libraries, the
> Qt libraries, and the compiler and platform runtime libraries listed in
> [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md), and to convey the resulting
> work. You must comply with the GNU General Public License in all respects for all of
> the code used other than those libraries. If you modify this file, you may extend
> this exception to your version of the file, but you are not obliged to do so.

## Why this exception exists

The pre-built Windows installer is produced with PyInstaller, which packages the
application together with every library it imports. Those libraries arrive as
pre-compiled wheels from the Python Package Index, and two of them carry terms the
Free Software Foundation regards as incompatible with the GPL:

- **Ipopt** (`ipopt-3.dll`), shipped inside the PySCIPOpt wheel, is distributed under
  the **Eclipse Public License**.
- **Intel compiler and OpenMP runtimes** (`libiomp5md.dll`, `libifcoremd.dll`,
  `libmmd.dll`, `svml_dispmd.dll`), also inside that wheel, are **proprietary
  redistributables**.

At the same time the application links GLPK, through `swiglpk`, which is GPLv3. A
single binary containing all of these is in tension no matter which licence the
application itself chooses. The conflict lies between the third-party components, not
with GSM ToolBox.

Section 7 of the GPL exists for exactly this case. The exception above removes the
ambiguity for the pre-built binary while leaving every other GPL obligation intact:
the source remains available, modifications remain copyleft, and the patent grant and
anti-tivoisation terms are unchanged.

## If you would rather avoid the question entirely

Install from source (see [`INSTALL.md`](INSTALL.md)). Nothing in this repository is
covered by the Eclipse Public License or by any proprietary licence; the third-party
libraries are then fetched by `pip` onto your own machine, which is not redistribution
and places no obligation on you.

You can also avoid the affected component specifically. The Ipopt and Intel libraries
arrive with PySCIPOpt, which provides the SCIP solver used for mixed-integer strain
design. Omit `pyscipopt` from the install and the application falls back to GLPK for
those analyses. That is slower on large problems and fully functional.
