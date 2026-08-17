"""Project / session state for GSM ToolBox.

A :class:`Project` wraps the working :class:`cobra.Model` together with:

* a pristine copy of the originally loaded model (for "diff vs original"),
* an undo/redo history (snapshot-based),
* attached datasets, cached analysis results and free-text notes,
* save/load to a single ``.gsmtbx`` file (a zip of the model + a JSON manifest).

The GUI never edits the model directly; it goes through :meth:`Project.apply_edit`
so that every change is snapshotted and reversible.
"""

from __future__ import annotations

import json
import os
import tempfile
import zipfile
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

import cobra

from . import io_models
from .categories import CategoryManager
from .flux_state import StrategyManager

# How many undo steps to retain. Snapshots are full model copies, so we cap this
# to keep memory bounded on genome-scale models.
_MAX_HISTORY = 25

_MANIFEST_NAME = "manifest.json"
_MODEL_NAME = "model.xml"
PROJECT_EXT = ".gsmtbx"


class ProjectError(Exception):
    """Raised on project save/load problems."""


@dataclass
class Project:
    """Holds the working model and all session state."""

    model: cobra.Model
    original_model: cobra.Model
    source_path: Optional[str] = None     # where the model was imported from
    project_path: Optional[str] = None    # the .gsmtbx file, once saved
    notes: str = ""
    datasets: Dict[str, Any] = field(default_factory=dict)
    results: Dict[str, Any] = field(default_factory=dict)
    settings: Dict[str, Any] = field(default_factory=dict)

    categories: CategoryManager = field(default_factory=CategoryManager)
    strategies: StrategyManager = field(default_factory=StrategyManager)

    _undo: List[cobra.Model] = field(default_factory=list, repr=False)
    _redo: List[cobra.Model] = field(default_factory=list, repr=False)
    _dirty: bool = False

    # ----- construction -------------------------------------------------
    @classmethod
    def from_model_file(cls, path: str) -> "Project":
        """Create a project by importing a model file (SBML/JSON/MAT)."""
        model = io_models.load_model(path)
        return cls(model=model, original_model=model.copy(), source_path=path)

    @classmethod
    def from_model(cls, model: cobra.Model, source_path: Optional[str] = None) -> "Project":
        return cls(model=model, original_model=model.copy(), source_path=source_path)

    # ----- editing with undo -------------------------------------------
    def apply_edit(self, edit_fn: Callable[[cobra.Model], Any]) -> Any:
        """Snapshot the model, then apply ``edit_fn(model)``.

        On any exception the snapshot is restored so the model is never left in a
        half-edited state. Returns whatever ``edit_fn`` returns.
        """
        snapshot = self.model.copy()
        try:
            result = edit_fn(self.model)
        except Exception:
            # Roll back to the snapshot we just took.
            self.model = snapshot
            raise
        self._undo.append(snapshot)
        if len(self._undo) > _MAX_HISTORY:
            self._undo.pop(0)
        self._redo.clear()
        self._dirty = True
        return result

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo(self) -> None:
        if not self._undo:
            return
        self._redo.append(self.model.copy())
        self.model = self._undo.pop()
        self._dirty = True

    def redo(self) -> None:
        if not self._redo:
            return
        self._undo.append(self.model.copy())
        self.model = self._redo.pop()
        self._dirty = True

    # ----- diff vs original --------------------------------------------
    def diff(self) -> Dict[str, list]:
        """Compare the working model to the original; report structural changes."""
        orig_rxns = {r.id for r in self.original_model.reactions}
        cur_rxns = {r.id for r in self.model.reactions}
        added = sorted(cur_rxns - orig_rxns)
        removed = sorted(orig_rxns - cur_rxns)

        changed_bounds = []
        for rid in sorted(cur_rxns & orig_rxns):
            o = self.original_model.reactions.get_by_id(rid)
            c = self.model.reactions.get_by_id(rid)
            if (o.lower_bound, o.upper_bound) != (c.lower_bound, c.upper_bound):
                changed_bounds.append(rid)
        orig_mets = {m.id for m in self.original_model.metabolites}
        cur_mets = {m.id for m in self.model.metabolites}
        return {
            "added_reactions": added,
            "removed_reactions": removed,
            "changed_bounds": changed_bounds,
            "added_metabolites": sorted(cur_mets - orig_mets),
            "removed_metabolites": sorted(orig_mets - cur_mets),
        }

    @property
    def is_modified(self) -> bool:
        return self._dirty

    def mark_saved(self) -> None:
        self._dirty = False

    # ----- persistence (.gsmtbx) ---------------------------------------
    def save(self, path: Optional[str] = None) -> str:
        """Save the project to a ``.gsmtbx`` file. Returns the path written."""
        path = path or self.project_path
        if not path:
            raise ProjectError("No project path provided.")
        if not path.lower().endswith(PROJECT_EXT):
            path += PROJECT_EXT

        manifest = {
            "format": "gsmtbx",
            "version": 1,
            "source_path": self.source_path,
            "notes": self.notes,
            "settings": self.settings,
            "categories": self.categories.to_list(),
            "strategies": self.strategies.to_list(),
            "diff": self.diff(),
        }
        try:
            with tempfile.TemporaryDirectory() as tmp:
                model_tmp = os.path.join(tmp, _MODEL_NAME)
                io_models.save_model(self.model, model_tmp)
                with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
                    zf.write(model_tmp, _MODEL_NAME)
                    zf.writestr(_MANIFEST_NAME, json.dumps(manifest, indent=2))
        except Exception as exc:  # noqa: BLE001
            raise ProjectError(f"Could not save project:\n{exc}") from exc

        self.project_path = path
        self.mark_saved()
        return path

    @classmethod
    def load(cls, path: str) -> "Project":
        """Load a project from a ``.gsmtbx`` file."""
        if not os.path.exists(path):
            raise ProjectError(f"Project file not found:\n{path}")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                with zipfile.ZipFile(path, "r") as zf:
                    zf.extractall(tmp)
                model_path = os.path.join(tmp, _MODEL_NAME)
                model = io_models.load_model(model_path)
                manifest_path = os.path.join(tmp, _MANIFEST_NAME)
                manifest = {}
                if os.path.exists(manifest_path):
                    with open(manifest_path, "r", encoding="utf-8") as fh:
                        manifest = json.load(fh)
        except ProjectError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise ProjectError(f"Could not open project:\n{exc}") from exc

        project = cls(
            model=model,
            original_model=model.copy(),
            source_path=manifest.get("source_path"),
            project_path=path,
            notes=manifest.get("notes", ""),
            settings=manifest.get("settings", {}),
        )
        project.categories.load_list(manifest.get("categories", []))
        project.strategies.load_list(manifest.get("strategies", []))
        return project
