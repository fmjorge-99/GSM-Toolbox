"""Persistent on-disk cache so the app stays useful offline.

Everything the app downloads is stored under ``~/.gsm_toolbox`` and recalled on later
runs, with no internet needed. That covers reaction databases, the MetaNetX
cross-reference files and molecule structure images. The Settings tab lets the user
inspect and clear the cache by category.

**Portable mode.** If a directory named ``.gsm_toolbox`` sits beside the executable, it
is used instead of the one in the user profile. That is what makes the portable Windows
build genuinely portable: the databases travel with it on a USB stick, and a machine the
user cannot install software on leaves nothing behind. The environment variable
``GSM_TOOLBOX_HOME`` overrides both, for a shared or scripted install.
"""

from __future__ import annotations

import hashlib
import os
import sys
from typing import Dict, List


def _app_dir() -> str:
    """The folder the application was launched from.

    Under PyInstaller ``sys.executable`` is the frozen exe, which is what a portable
    copy sits next to. From source it is the interpreter, so the package location is
    used instead and a stray ``.gsm_toolbox`` in a virtual environment cannot be picked
    up by accident.
    """
    if getattr(sys, "frozen", False):
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _resolve_base() -> str:
    override = os.environ.get("GSM_TOOLBOX_HOME", "").strip()
    if override:
        return os.path.abspath(os.path.expanduser(override))
    # Opt-in, and only when the folder already exists. Creating it automatically would
    # scatter caches beside every checkout and quietly strand a user's downloads.
    portable = os.path.join(_app_dir(), ".gsm_toolbox")
    if os.path.isdir(portable):
        return portable
    return os.path.join(os.path.expanduser("~"), ".gsm_toolbox")


_BASE = _resolve_base()


def base_dir() -> str:
    os.makedirs(_BASE, exist_ok=True)
    return _BASE


def is_portable() -> bool:
    """True when the cache lives beside the application rather than in the profile."""
    return os.path.normcase(_BASE) != os.path.normcase(
        os.path.join(os.path.expanduser("~"), ".gsm_toolbox"))


def _sub(name: str) -> str:
    path = os.path.join(base_dir(), name)
    os.makedirs(path, exist_ok=True)
    return path


def databases_dir() -> str:
    return _sub("databases")


def structures_dir() -> str:
    return _sub("structures")


# ---- molecule-structure image cache --------------------------------------
def _key_hash(key: str) -> str:
    return hashlib.sha1(key.encode("utf-8")).hexdigest()


def image_path(key: str) -> str:
    return os.path.join(structures_dir(), f"{_key_hash(key)}.png")


def get_image(key: str) -> bytes | None:
    """Return cached PNG bytes for ``key`` (a stable descriptor), or None."""
    path = image_path(key)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            with open(path, "rb") as fh:
                return fh.read()
        except OSError:
            return None
    return None


def put_image(key: str, data: bytes) -> None:
    if not data:
        return
    try:
        with open(image_path(key), "wb") as fh:
            fh.write(data)
    except OSError:
        pass


# ---- inspection / clearing (for the Settings tab) ------------------------
# Category label -> (directory, description)
CATEGORIES: Dict[str, str] = {
    "Reaction databases": "databases",
    "Molecule structures": "structures",
}


def _category_dir(category: str) -> str:
    return _sub(CATEGORIES.get(category, category))


def list_category(category: str) -> List[dict]:
    """List cached files in a category: ``[{path, name, size}]`` (largest first)."""
    d = _category_dir(category)
    out = []
    for name in os.listdir(d) if os.path.isdir(d) else []:
        path = os.path.join(d, name)
        if os.path.isfile(path):
            try:
                out.append({"path": path, "name": name, "size": os.path.getsize(path)})
            except OSError:
                continue
    out.sort(key=lambda r: r["size"], reverse=True)
    return out


def category_total(category: str) -> int:
    return sum(f["size"] for f in list_category(category))


def delete_paths(paths: List[str]) -> int:
    """Delete the given files; return the number of bytes freed."""
    freed = 0
    for p in paths:
        try:
            if os.path.isfile(p):
                freed += os.path.getsize(p)
                os.remove(p)
        except OSError:
            continue
    return freed


def clear_category(category: str) -> int:
    return delete_paths([f["path"] for f in list_category(category)])


def human_size(n: int) -> str:
    step = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if step < 1024 or unit == "TB":
            return f"{step:.0f} {unit}" if unit == "B" else f"{step:.1f} {unit}"
        step /= 1024
    return f"{n} B"
