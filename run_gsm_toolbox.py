"""Frozen-app entry point.

PyInstaller runs its entry script as ``__main__`` with no package context, which
breaks the relative imports inside ``gsm_toolbox/app.py``. This thin launcher
imports the package properly so those relative imports resolve.

Crucially, it calls ``multiprocessing.freeze_support()`` BEFORE importing the GUI.
Libraries such as StrainDesign use multiprocessing; on Windows (spawn start
method) each worker re-launches the frozen exe. Without freeze_support() those
workers would re-run the whole application and open extra windows. freeze_support()
makes a spawned process behave as a worker and exit, so only one window appears.
"""

import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from gsm_toolbox.app import main

    sys.exit(main())
