"""Loading and saving genome-scale metabolic models.

Supports the formats COBRApy understands: SBML (``.xml``/``.sbml``), JSON
(``.json``) and MATLAB (``.mat``). All loaders raise :class:`ModelLoadError`
with a human-readable message so the GUI can show a friendly dialog instead of
a raw traceback.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import List

import cobra
from cobra.io import (
    load_json_model,
    load_matlab_model,
    model_from_dict,
    read_sbml_model,
    save_json_model,
    save_matlab_model,
    write_sbml_model,
)

# Map of lower-cased file extension -> (loader, saver). Savers may be ``None``
# when we only support reading that format.
_SBML_EXTS = {".xml", ".sbml"}
_JSON_EXTS = {".json"}
_MAT_EXTS = {".mat"}

SUPPORTED_LOAD_EXTS = _SBML_EXTS | _JSON_EXTS | _MAT_EXTS
SUPPORTED_SAVE_EXTS = _SBML_EXTS | _JSON_EXTS | _MAT_EXTS


class ModelLoadError(Exception):
    """Raised when a model file cannot be read or parsed."""


class ModelSaveError(Exception):
    """Raised when a model cannot be written to disk."""


@dataclass
class ModelSummary:
    """A quick, human-readable snapshot of a model's contents."""

    model_id: str
    name: str
    n_reactions: int
    n_metabolites: int
    n_genes: int
    n_exchanges: int
    compartments: dict = field(default_factory=dict)
    objective: str = ""

    def as_rows(self) -> List[tuple]:
        """Return ``(label, value)`` rows for display in a summary table."""
        return [
            ("Model ID", self.model_id or "—"),
            ("Name", self.name or "—"),
            ("Reactions", self.n_reactions),
            ("Metabolites", self.n_metabolites),
            ("Genes", self.n_genes),
            ("Exchange reactions", self.n_exchanges),
            ("Compartments", ", ".join(f"{k}: {v}" for k, v in self.compartments.items()) or "—"),
            ("Objective", self.objective or "—"),
        ]


def _ext(path: str) -> str:
    return os.path.splitext(path)[1].lower()


def _coerce_annotation(ann):
    """Normalise a cobra-JSON ``annotation`` field into a plain dict.

    Some databases (notably the BiGG *universal* model) store annotations as a
    list of ``[key, value]`` pairs instead of the ``{key: value}`` mapping that
    cobra requires. Values for a repeated key are collected into a list.
    """
    if isinstance(ann, dict) or ann is None:
        return ann
    if isinstance(ann, (list, tuple)):
        out: dict = {}
        for item in ann:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                key, value = item
                if key in out:
                    existing = out[key]
                    if isinstance(existing, list):
                        existing.append(value)
                    else:
                        out[key] = [existing, value]
                else:
                    out[key] = value
        return out
    # Anything else (str, etc.) -> drop it, cobra would reject it.
    return {}


def _sanitize_cobra_dict(data: dict) -> dict:
    """Make a raw cobra-JSON ``dict`` safe for :func:`model_from_dict`.

    Fixes list-style annotations and missing compartment metadata so that
    BiGG-format universal models load cleanly.
    """
    for section in ("metabolites", "reactions", "genes"):
        for entry in data.get(section, []) or []:
            if not isinstance(entry, dict):
                continue
            if "annotation" in entry:
                entry["annotation"] = _coerce_annotation(entry["annotation"])

    # cobra needs every metabolite to name a compartment, and a top-level
    # compartments mapping. The BiGG universal model leaves both empty.
    compartments = data.get("compartments")
    if not isinstance(compartments, dict):
        compartments = {}
        data["compartments"] = compartments
    for met in data.get("metabolites", []) or []:
        if not isinstance(met, dict):
            continue
        comp = met.get("compartment")
        if not comp:
            # Infer from a trailing "_<comp>" suffix on the id (BiGG style).
            mid = str(met.get("id", ""))
            comp = mid.rsplit("_", 1)[-1] if "_" in mid else ""
            if not comp:
                comp = "c"
            met["compartment"] = comp
        compartments.setdefault(comp, comp)
    return data


def _load_json_tolerant(path: str) -> cobra.Model:
    """Load a JSON model, repairing BiGG-style schema quirks on failure."""
    try:
        return load_json_model(path)
    except (TypeError, KeyError, ValueError):
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return model_from_dict(_sanitize_cobra_dict(data))


def load_model(path: str) -> cobra.Model:
    """Load a metabolic model from ``path``, dispatching on file extension.

    Raises :class:`ModelLoadError` with a friendly message on any failure.
    """
    if not os.path.exists(path):
        raise ModelLoadError(f"File not found:\n{path}")

    ext = _ext(path)
    try:
        if ext in _SBML_EXTS:
            model = read_sbml_model(path)
        elif ext in _JSON_EXTS:
            model = _load_json_tolerant(path)
        elif ext in _MAT_EXTS:
            model = load_matlab_model(path)
        else:
            raise ModelLoadError(
                f"Unsupported file type '{ext}'.\n"
                "Please choose an SBML (.xml/.sbml), JSON (.json) or MATLAB (.mat) model."
            )
    except ModelLoadError:
        raise
    except Exception as exc:  # noqa: BLE001 - we deliberately wrap everything
        raise ModelLoadError(
            f"Could not read the model file:\n{os.path.basename(path)}\n\nDetails: {exc}"
        ) from exc

    if not model.reactions:
        raise ModelLoadError(
            "The file was read but contains no reactions. "
            "It may not be a valid metabolic model."
        )
    return model


def save_model(model: cobra.Model, path: str) -> None:
    """Save ``model`` to ``path``, dispatching on file extension."""
    ext = _ext(path)
    try:
        if ext in _SBML_EXTS:
            write_sbml_model(model, path)
        elif ext in _JSON_EXTS:
            save_json_model(model, path)
        elif ext in _MAT_EXTS:
            save_matlab_model(model, path)
        else:
            raise ModelSaveError(
                f"Unsupported file type '{ext}'. Use .xml, .sbml, .json or .mat."
            )
    except ModelSaveError:
        raise
    except Exception as exc:  # noqa: BLE001
        raise ModelSaveError(f"Could not save the model:\n{exc}") from exc


def summarize(model: cobra.Model) -> ModelSummary:
    """Build a :class:`ModelSummary` for ``model``."""
    try:
        objective = str(model.objective.expression)
    except Exception:  # noqa: BLE001
        objective = ""
    return ModelSummary(
        model_id=model.id or "",
        name=getattr(model, "name", "") or "",
        n_reactions=len(model.reactions),
        n_metabolites=len(model.metabolites),
        n_genes=len(model.genes),
        n_exchanges=len(model.exchanges),
        compartments=dict(model.compartments),
        objective=objective,
    )
