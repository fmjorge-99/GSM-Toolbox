"""Allow ``python -m gsm_toolbox`` to launch the app."""

import sys

from .app import main

if __name__ == "__main__":
    sys.exit(main())
