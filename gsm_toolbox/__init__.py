"""GSM ToolBox — a user-friendly desktop app for genome-scale metabolic modeling.

Copyright (C) 2026 Jorge Fernández Méndez.

This program is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. It is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details. You should have received a copy of the licence along
with this program; if not, see <https://www.gnu.org/licenses/>.

An additional permission under section 7, covering the solver and platform runtime
libraries bundled in the pre-built installer, is stated in LICENSE-EXCEPTION.md.

The package is split into two layers:

* :mod:`gsm_toolbox.core` — a pure-Python science engine that wraps COBRApy and
  related libraries. It contains **no Qt imports** and is fully testable headlessly.
* :mod:`gsm_toolbox.gui` — the PySide6 desktop interface that drives the core engine.
"""

__version__ = "0.3.10"
__app_name__ = "GSM ToolBox"
