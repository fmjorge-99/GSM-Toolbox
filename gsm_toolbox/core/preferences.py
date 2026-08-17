"""User preferences, persisted as JSON in ``~/.gsm_toolbox/preferences.json``.

Deliberately tiny: a flat key/value store read on demand and written atomically. It holds
the switches that decide which optional parts of the app are visible at all — most
importantly whether the thermodynamics (MDF) suite is enabled, since that feature needs a
large external dataset and is hidden until the user opts in.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict

from . import cache

# --- keys -----------------------------------------------------------------------------
MDF_ENABLED = "mdf_enabled"            # bool — show the thermodynamics (MDF) suite
RETRORULES_SEED = "retrorules_seed"    # int  — deterministic seed for rule search
RETRORULES_SEED_DEFAULT = 12345
SELENZYME_ENABLED = "selenzyme_enabled"   # bool — reaction-similarity enzyme search

# --- General --------------------------------------------------------------------------
RESULTS_DIR = "results_dir"                  # str  — auto-save analysis tables here
AUTOSAVE_RESULTS = "autosave_results"        # bool — write a CSV per analysis
CONFIRM_ON_EXIT = "confirm_on_exit"          # bool — ask before closing with unsaved work
RESTORE_LAST_PROJECT = "restore_last_project"

# --- Appearance -----------------------------------------------------------------------
FONT_SIZE = "font_size"                      # int  — base point size
INFO_PANEL_POSITION = "info_panel_position"  # "right" | "bottom"
TABLE_DENSITY = "table_density"              # "comfortable" | "compact"
SHOW_EXPLORER = "show_explorer"
SHOW_CATEGORIES = "show_categories"
SHOW_INFO = "show_info"

# --- Pathway Finder -------------------------------------------------------------------
RR_DIAMETER = "retrorules_diameter"          # int  — rule specificity
RR_PLAUSIBILITY = "retrorules_plausibility"  # bool — hide chemically implausible routes
STRUCTURE_FETCH = "structure_fetch_online"   # bool — resolve structures from the web
SEARCH_ALGORITHM = "search_algorithm"        # "retro" | "expansion"
SEARCH_MAX_STEPS = "search_max_steps"
SEARCH_ALTERNATIVES = "search_alternatives"
FLUX_CARRYING_STARTS = "flux_carrying_starts"

# --- Analysis & solver ----------------------------------------------------------------
SOLVER = "solver"                            # "" = let COBRApy choose
FVA_FRACTION = "fva_fraction"
WORKER_PROCESSES = "worker_processes"
PRODUCTION_GROWTH_FLOOR = "production_growth_floor"

# --- Regulation & dynamics ------------------------------------------------------------
REGULATION_ENABLED = "regulation_enabled"
REGULATION_RULESET = "regulation_ruleset"    # str — path to a rule-set JSON ("" = bundled)
REGULATION_ACTIVATION = "regulation_activation_threshold"
DFBA_STEP_H = "dfba_step_h"
DFBA_DURATION_H = "dfba_duration_h"

# --- Data & storage -------------------------------------------------------------------
OFFLINE_MODE = "offline_mode"                # bool — never reach the network
ALLOW_DOWNLOADS = "allow_downloads"          # bool — consent for one-off reference data

# --- Shortcuts ------------------------------------------------------------------------
TOOLBAR_ACTIONS = "toolbar_actions"          # list[str] — ids shown on the quick-access bar
#: The default quick-access bar. "Growth Settings" and "Network Visualization" are here
#: rather than as main tabs, which is what freed the tab strip for Dynamic Analysis.
TOOLBAR_DEFAULT = ["open_model", "save_project", "add_reaction", "undo", "redo",
                   "growth_settings", "network_visualization"]

_DEFAULTS: Dict[str, Any] = {
    MDF_ENABLED: False,                # opt-in: hidden until the user enables it
    RETRORULES_SEED: RETRORULES_SEED_DEFAULT,
    SELENZYME_ENABLED: False,          # opt-in: needs a one-off reference-data download

    RESULTS_DIR: "",
    AUTOSAVE_RESULTS: False,
    CONFIRM_ON_EXIT: True,
    RESTORE_LAST_PROJECT: False,

    FONT_SIZE: 9,
    INFO_PANEL_POSITION: "right",
    TABLE_DENSITY: "comfortable",
    SHOW_EXPLORER: True,
    SHOW_CATEGORIES: True,
    #: Off by default. The Information dock is a detail view for whatever is selected,
    #: so it is empty on startup while taking a column of width from the tables and maps
    #: that are not. It fills the moment anything is clicked, and View ▸ Information
    #: turns it on for good.
    SHOW_INFO: False,

    RR_DIAMETER: 6,
    RR_PLAUSIBILITY: True,
    STRUCTURE_FETCH: True,
    SEARCH_ALGORITHM: "retro",
    SEARCH_MAX_STEPS: 25,
    SEARCH_ALTERNATIVES: 1,
    FLUX_CARRYING_STARTS: False,       # off by default: preserves earlier results

    SOLVER: "",
    FVA_FRACTION: 0.9,
    WORKER_PROCESSES: 1,               # >1 is unsafe with spawn on Windows
    PRODUCTION_GROWTH_FLOOR: 0.1,

    REGULATION_ENABLED: False,         # opt-in until the rule thresholds are sourced
    REGULATION_RULESET: "",
    REGULATION_ACTIVATION: 0.05,
    DFBA_STEP_H: 4.0,
    DFBA_DURATION_H: 96.0,

    OFFLINE_MODE: False,
    ALLOW_DOWNLOADS: True,

    TOOLBAR_ACTIONS: list(TOOLBAR_DEFAULT),
}

# --- migrations -----------------------------------------------------------------------
#: Bumped when a stored preference needs correcting rather than merely defaulting
#: differently. A changed default only reaches a fresh install; anyone who has opened
#: Preferences once has the old value written to disk.
SCHEMA_VERSION = "preferences_version"
_CURRENT_SCHEMA = 2


def _migrate(data: Dict[str, Any]) -> bool:
    """Bring a stored preferences file up to date. Returns True if anything changed."""
    version = int(data.get(SCHEMA_VERSION, 1) or 1)
    changed = False

    if version < 2:
        # The Information dock used to default to visible. It is a detail view for
        # whatever is selected, so at startup it is empty while taking a column of width
        # from the tables and maps that are not — and the stored `True` came from that
        # old default rather than from anyone choosing it. Re-enable from
        # View ▸ Information; the choice then persists.
        if data.get(SHOW_INFO) is True:
            data[SHOW_INFO] = False
            changed = True

    if version != _CURRENT_SCHEMA:
        data[SCHEMA_VERSION] = _CURRENT_SCHEMA
        changed = True
    return changed

_CACHE: Dict[str, Any] | None = None


def _path() -> str:
    return os.path.join(cache.base_dir(), "preferences.json")


def _load() -> Dict[str, Any]:
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    data = dict(_DEFAULTS)
    try:
        with open(_path()) as fh:
            stored = json.load(fh)
        if isinstance(stored, dict):
            data.update(stored)
    except Exception:  # noqa: BLE001 — a missing/corrupt file just means defaults
        pass
    if _migrate(data):
        _CACHE = data
        _write(data)
    _CACHE = data
    return data


def get(key: str, default: Any = None) -> Any:
    data = _load()
    return data.get(key, _DEFAULTS.get(key, default))


def _write(data: Dict[str, Any]) -> None:
    try:
        os.makedirs(cache.base_dir(), exist_ok=True)
        tmp = _path() + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(data, fh, indent=1)
        os.replace(tmp, _path())        # atomic: never leave a half-written file
    except Exception:  # noqa: BLE001 — preferences must never break the app
        pass


def set(key: str, value: Any) -> None:  # noqa: A001 — mirrors dict/QSettings naming
    data = _load()
    data[key] = value
    _write(data)


def reload() -> None:
    """Drop the in-memory copy so the next read re-reads the file."""
    global _CACHE
    _CACHE = None


# --- convenience ----------------------------------------------------------------------
def mdf_enabled() -> bool:
    return bool(get(MDF_ENABLED, False))


def selenzyme_enabled() -> bool:
    return bool(get(SELENZYME_ENABLED, False))


def retrorules_seed() -> int:
    try:
        return int(get(RETRORULES_SEED, RETRORULES_SEED_DEFAULT))
    except (TypeError, ValueError):
        return RETRORULES_SEED_DEFAULT
