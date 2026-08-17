"""Guarantee that ``sys.stdout`` and ``sys.stderr`` are writable.

A PyInstaller **windowed** build (``console=False``) starts with ``sys.stdout`` and
``sys.stderr`` set to ``None`` — there is no console to write to. Python's own ``print``
and the ``logging`` module cope with that, but third-party libraries generally do not.
``tqdm`` in particular writes to ``sys.stderr`` unconditionally, so the moment any
dependency draws a progress bar the frozen app dies with:

    'NoneType' object has no attribute 'write'

That is exactly the error the thermodynamics analysis produced: eQuilibrator draws a
tqdm bar while decomposing compounds it cannot resolve by accession — the novel
chemistry of a cannabinoid route — so the crash appeared only for the routes where the
analysis was most needed, and only in the packaged app, never in a development run.

Rather than silencing the output, this redirects it to a rotating log file in the user's
cache directory. A packaged app that swallows its diagnostics is much harder to support,
and the file is the only place a traceback can go once there is no console.
"""

from __future__ import annotations

import os
import sys
from typing import Optional

#: Kept module-level so the handle stays open for the process lifetime.
_LOG_HANDLE = None
_LOG_PATH: Optional[str] = None

#: Roll the log over rather than letting it grow without bound.
_MAX_BYTES = 2 * 1024 * 1024


class _NullStream:
    """Absorbs writes when even the log file cannot be opened (read-only install)."""

    def write(self, _data):          # noqa: D401 - stream protocol
        return 0

    def flush(self):
        return None

    def isatty(self):
        return False

    def fileno(self):
        raise OSError("no file descriptor")

    def close(self):
        return None


def log_path() -> str:
    """Where console output is redirected when there is no console."""
    from . import cache
    return os.path.join(cache.base_dir(), "logs", "console.log")


def _open_log():
    path = log_path()
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path) and os.path.getsize(path) > _MAX_BYTES:
            backup = path + ".1"
            try:
                if os.path.exists(backup):
                    os.remove(backup)
                os.replace(path, backup)
            except OSError:
                pass                 # rotation is a nicety, never a reason to fail
        handle = open(path, "a", encoding="utf-8", errors="replace", buffering=1)
    except Exception:  # noqa: BLE001 - read-only install, locked file, anything
        return None, None
    return handle, path


def install() -> Optional[str]:
    """Replace any ``None`` standard stream. Returns the log path, or None.

    Safe to call more than once and safe to call when the streams are already fine —
    a normal ``python -m gsm_toolbox`` run is left completely untouched.
    """
    global _LOG_HANDLE, _LOG_PATH
    if sys.stdout is not None and sys.stderr is not None:
        return None                  # running with a console: nothing to do

    if _LOG_HANDLE is None:
        _LOG_HANDLE, _LOG_PATH = _open_log()
    stream = _LOG_HANDLE or _NullStream()

    if sys.stdout is None:
        sys.stdout = stream
    if sys.stderr is None:
        sys.stderr = stream
    # __stdout__/__stderr__ are None too in a windowed build, and some libraries reach
    # for them directly to bypass a redirected sys.stdout.
    if getattr(sys, "__stdout__", None) is None:
        sys.__stdout__ = stream
    if getattr(sys, "__stderr__", None) is None:
        sys.__stderr__ = stream
    return _LOG_PATH
