"""Let data pickled by numpy 2 load under numpy 1, and vice versa.

numpy 2 renamed the internal package ``numpy.core`` to ``numpy._core``. Anything pickled
by numpy 2 — including the ``.npz`` of component-contribution parameters that eQuilibrator
downloads — records the *new* path, and unpickling it under numpy 1 fails with::

    No module named 'numpy._core'

numpy 1.26 anticipates this and ships a ``numpy/_core/`` shim for exactly this purpose.
But nothing imports it: pickle reaches it dynamically through ``find_class``, so a
static-analysis bundler (PyInstaller) does not see it and leaves it out of the frozen
app. The result is the worst kind of bug — the thermodynamics analysis works perfectly in
development and dies in the packaged app, after a long download, with an error that names
numpy rather than anything the user did.

The spec now bundles ``numpy`` wholesale, which is the real fix. This module is the safety
net: it maps any missing ``numpy._core*`` name onto its ``numpy.core*`` equivalent (and
the reverse, for a numpy-1 pickle read under numpy 2), so an incomplete build degrades to
"slightly slower import" rather than a hard failure.
"""

from __future__ import annotations

import importlib
import sys
from typing import List

#: The submodules pickled numpy data actually reaches for.
_SUBMODULES = ("multiarray", "umath", "_multiarray_umath", "_internal", "_dtype",
               "_dtype_ctypes", "numeric", "_methods")


def install() -> List[str]:
    """Alias the missing numpy internal package onto the one this build has.

    Returns the alias names that had to be created — empty when numpy is complete, which
    is the normal case. Safe to call repeatedly and safe when numpy is absent.
    """
    try:
        import numpy
    except Exception:  # noqa: BLE001 - no numpy: nothing to reconcile
        return []

    created: List[str] = []
    # numpy 2 layout is "_core"; numpy 1 layout is "core". Whichever this build has is
    # the source; the other name is the alias we may need to provide.
    for missing, present in (("numpy._core", "numpy.core"),
                             ("numpy.core", "numpy._core")):
        if _importable(missing) or not _importable(present):
            continue
        source = importlib.import_module(present)
        sys.modules[missing] = source
        created.append(missing)
        for sub in _SUBMODULES:
            try:
                mod = importlib.import_module(f"{present}.{sub}")
            except Exception:  # noqa: BLE001 - not every submodule exists in every version
                continue
            sys.modules[f"{missing}.{sub}"] = mod
            created.append(f"{missing}.{sub}")
    return created


def _importable(name: str) -> bool:
    if name in sys.modules:
        return True
    try:
        importlib.import_module(name)
        return True
    except Exception:  # noqa: BLE001 - ImportError and anything a broken build raises
        return False
