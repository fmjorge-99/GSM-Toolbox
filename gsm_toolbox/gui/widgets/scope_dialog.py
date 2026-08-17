"""A small reusable multi-select dialog for choosing one or more categories or
subsystems to focus a map on. Shared by the Network Map and the Escher Visualizer
so both offer the identical Focus ▸ Category / Subsystem experience."""

from __future__ import annotations

from typing import Dict, List, Optional, Set

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
)


def choose_scope(parent, kind: str, items: Dict[str, List[str]],
                 selected: Set[str]) -> Optional[Set[str]]:
    """Show a checkable list of ``items`` (name -> reaction ids). ``kind`` is the noun
    ("Category"/"Subsystem"). Returns the new set of ticked names, or ``None`` if the
    user cancelled."""
    noun = kind.lower()
    dlg = QDialog(parent)
    dlg.setWindowTitle(f"Choose {noun}(s) to map")
    lay = QVBoxLayout(dlg)
    lay.addWidget(QLabel(f"Tick one or more {noun}s to draw. Their reactions are combined "
                         "into a single focused map."))
    lst = QListWidget()
    lst.setSelectionMode(QAbstractItemView.NoSelection)
    for name in sorted(items):
        it = QListWidgetItem(f"{name}  ({len(items[name])})")
        it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
        it.setCheckState(Qt.Checked if name in selected else Qt.Unchecked)
        it.setData(Qt.UserRole, name)
        lst.addItem(it)
    lay.addWidget(lst, 1)
    bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    bb.accepted.connect(dlg.accept)
    bb.rejected.connect(dlg.reject)
    lay.addWidget(bb)
    dlg.resize(360, 460)
    if dlg.exec() != QDialog.Accepted:
        return None
    return {lst.item(i).data(Qt.UserRole) for i in range(lst.count())
            if lst.item(i).checkState() == Qt.Checked}


def model_subsystems(model) -> Dict[str, List[str]]:
    """{subsystem name: [reaction ids]} for a model (skips reactions with no subsystem)."""
    out: Dict[str, List[str]] = {}
    if model is None:
        return out
    for r in model.reactions:
        sub = (getattr(r, "subsystem", "") or "").strip()
        if sub:
            out.setdefault(sub, []).append(r.id)
    return out
