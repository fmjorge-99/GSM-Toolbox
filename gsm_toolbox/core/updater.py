"""Check (never auto-install) the versions of the toolbox, its dependencies, and its
downloadable data, and tell the user what — if anything — is out of date and how to
update it by hand.

Deliberately read-only. A frozen desktop app cannot reliably pip-install into itself, and
silently mutating a scientist's environment is worse than telling them a one-line command.
So every check returns a *status* and a *suggested command*; acting on it is the user's
call. All network access is optional and short-timeout: offline simply yields "unknown".
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List, Optional

from . import cache

# Dependencies worth watching: the scientific/UI core plus the optional thermo add-on.
_TRACKED_DEPS = [
    "cobra", "PySide6", "rdkit", "numpy", "scipy", "pandas", "optlang",
    "python-libsbml", "escher", "equilibrator-api", "equilibrator-cache",
    "component-contribution", "requests", "matplotlib",
]

OK = "ok"
OUTDATED = "outdated"
MISSING = "missing"
UNKNOWN = "unknown"       # e.g. no network to check the latest version
INFO = "info"             # reported for context, no update concept


@dataclass
class Component:
    name: str
    kind: str                       # "dependency" | "data" | "app"
    current: str = ""
    latest: str = ""
    status: str = UNKNOWN
    detail: str = ""
    suggestion: str = ""            # a manual command or action, if any


@dataclass
class UpdateReport:
    components: List[Component] = field(default_factory=list)
    checked_online: bool = False

    @property
    def n_outdated(self) -> int:
        return sum(1 for c in self.components if c.status == OUTDATED)

    @property
    def n_missing(self) -> int:
        return sum(1 for c in self.components if c.status == MISSING)


def _installed_version(dist: str) -> Optional[str]:
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version(dist)
    except PackageNotFoundError:
        return None


def _pypi_latest(dist: str, *, timeout: float = 4.0) -> Optional[str]:
    """Latest release on PyPI, or None if offline / not found."""
    import json
    from urllib.request import urlopen
    try:
        with urlopen(f"https://pypi.org/pypi/{dist}/json", timeout=timeout) as fh:
            data = json.load(fh)
        return data.get("info", {}).get("version")
    except Exception:  # noqa: BLE001 — offline or PyPI hiccup: treat as unknown
        return None


def _cmp_versions(current: str, latest: str) -> int:
    """-1 if current < latest, 0 if equal, 1 if current > latest (best-effort)."""
    try:
        from packaging.version import Version
        a, b = Version(current), Version(latest)
        return -1 if a < b else (1 if a > b else 0)
    except Exception:  # noqa: BLE001 — fall back to a naive tuple compare
        def parts(v):
            return [int(x) for x in v.split(".") if x.isdigit()]
        pa, pb = parts(current), parts(latest)
        return -1 if pa < pb else (1 if pa > pb else 0)


def check_dependencies(*, online: bool = True) -> List[Component]:
    out: List[Component] = []
    for dist in _TRACKED_DEPS:
        cur = _installed_version(dist)
        c = Component(name=dist, kind="dependency", current=cur or "")
        if cur is None:
            # equilibrator-* are optional; the rest being absent is genuinely notable.
            optional = dist.startswith(("equilibrator", "component")) or dist == "escher"
            c.status = INFO if optional else MISSING
            c.detail = "not installed" + (" (optional)" if optional else "")
            c.suggestion = f"pip install {dist}"
            out.append(c)
            continue
        latest = _pypi_latest(dist) if online else None
        if latest is None:
            c.status = UNKNOWN
            c.detail = "installed; latest unknown (offline?)"
        else:
            c.latest = latest
            if _cmp_versions(cur, latest) < 0:
                c.status = OUTDATED
                c.suggestion = f"pip install -U {dist}"
            else:
                c.status = OK
        out.append(c)
    return out


def _dir_size(path: str) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


def _age_days(path: str) -> Optional[float]:
    import time
    try:
        return (time.time() - os.path.getmtime(path)) / 86400.0
    except OSError:
        return None


def check_data() -> List[Component]:
    """Downloadable data assets under ~/.gsm_toolbox: present? how big? how old?"""
    from .thermodynamics import equilibrator_available
    out: List[Component] = []
    base = cache.base_dir()

    # Reaction databases (BiGG / ModelSEED / merged / KEGG fetches).
    db_dir = cache.databases_dir()
    dbs = [f for f in os.listdir(db_dir)] if os.path.isdir(db_dir) else []
    json_dbs = [f for f in dbs if f.endswith(".json")]
    c = Component(name="Reaction databases", kind="data",
                  current=f"{len(json_dbs)} database file(s)")
    if json_dbs:
        c.status = INFO
        c.detail = f"{cache.human_size(_dir_size(db_dir))} in {db_dir}"
        c.suggestion = ("Refresh from the Pathway Design panel ▸ Fetch online, "
                        "or Manage databases ▸ re-merge.")
    else:
        c.status = MISSING
        c.detail = "no reaction database downloaded yet"
        c.suggestion = "Pathway Design ▸ Fetch online (BiGG universal / ModelSEED)."
    out.append(c)

    # RetroRules ruleset.
    rr_dir = os.path.join(base, "retrorules")
    rr_present = os.path.isdir(rr_dir) and bool(os.listdir(rr_dir))
    c = Component(name="RetroRules ruleset", kind="data")
    if rr_present:
        age = _age_days(rr_dir)
        c.current = "installed"
        c.status = INFO
        c.detail = (f"{cache.human_size(_dir_size(rr_dir))}"
                    + (f", ~{age:.0f} days old" if age is not None else ""))
        c.suggestion = ("Delete the ~/.gsm_toolbox/retrorules folder to force a fresh "
                        "download on next RetroRules run.")
    else:
        c.current = "not downloaded"
        c.status = INFO
        c.detail = "downloaded on demand the first time RetroRules Prediction is used"
    out.append(c)

    # eQuilibrator compound cache (thermodynamics).
    c = Component(name="eQuilibrator data (thermodynamics)", kind="data")
    if not equilibrator_available():
        c.current = "add-on not installed"
        c.status = INFO
        c.detail = "MDF thermodynamic analysis is unavailable without it"
        c.suggestion = "pip install equilibrator-api"
    else:
        cache_file = _equilibrator_cache_file()
        if cache_file:
            c.current = "downloaded"
            c.status = INFO
            try:
                c.detail = f"{cache.human_size(os.path.getsize(cache_file))} at {cache_file}"
            except OSError:
                c.detail = cache_file
        else:
            c.current = "installed; data not yet downloaded"
            c.status = INFO
            c.detail = "downloaded automatically the first time MDF analysis runs"
    out.append(c)
    return out


def _equilibrator_cache_file() -> Optional[str]:
    """Best-effort path to the downloaded eQuilibrator compound cache, if present."""
    candidates = [
        os.path.join(os.path.expanduser("~"), ".cache", "equilibrator"),
        os.path.join(os.environ.get("LOCALAPPDATA", ""), "equilibrator"),
    ]
    for d in candidates:
        if d and os.path.isdir(d):
            for root, _dirs, files in os.walk(d):
                for f in files:
                    if f.endswith(".sqlite"):
                        return os.path.join(root, f)
    return None


def check_app() -> Component:
    from .. import __version__
    c = Component(name="GSM ToolBox", kind="app", current=__version__, status=INFO)
    c.detail = "installed application version"
    c.suggestion = ("Updates ship as a new installer — reinstall the latest "
                    "GSM_ToolBox_Setup to upgrade the whole app at once.")
    return c


def build_report(*, online: bool = True) -> UpdateReport:
    """Assemble the full read-only update report."""
    r = UpdateReport(checked_online=online)
    r.components.append(check_app())
    deps = check_dependencies(online=online)
    r.components.extend(deps)
    r.components.extend(check_data())
    # We managed an online check iff at least one dependency resolved a latest version.
    r.checked_online = online and any(c.latest for c in deps)
    return r
