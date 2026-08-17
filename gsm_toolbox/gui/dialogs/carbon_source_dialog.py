"""Pick the organic carbon source(s) to feed, when switching to mixo- or heterotrophy.

Both modes are defined by having organic carbon available, but *which* carbon is a
decision the preset cannot make: a model may offer thirty candidate exchanges and the
answer depends on the experiment. Previously the preset reported what it had done and the
user then had to find the right exchange in the medium editor by hand.

This asks the question at the moment it arises, with the plausible answers already listed
and the uptake rate editable. Anything already open in the medium is pre-ticked, so
accepting without changes keeps the model as it is.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QHBoxLayout, QLabel,
    QLineEdit, QScrollArea, QVBoxLayout, QWidget)

from .. import style

#: Offered first: the substrates a cyanobacterial or general model is usually fed.
_COMMON = ("glc__D", "glyc", "ac", "succ", "fru", "xyl__D", "lac__L", "pyr", "cit",
           "mal__L", "akg", "for", "etoh", "arab__L", "gam", "man", "sucr", "malt")

#: Default uptake, mmol gDW⁻¹ h⁻¹. Nogales *et al.* use 0.85 for glucose in dark
#: heterotrophy and 0.38 mixotrophically; 0.85 is the safer starting point because it is
#: the rate the reference model was characterised at.
DEFAULT_UPTAKE = 0.85


class CarbonSourceDialog(QDialog):
    """Tick one or more organic carbon exchanges and set their uptake rate."""

    def __init__(self, parent, model, mode: str, notes: Sequence[str] = ()):
        super().__init__(parent)
        self.setWindowTitle(f"{mode.capitalize()} — choose a carbon source")
        self.resize(560, 560)
        from ..widgets.dialog_util import clamp_to_screen

        self._model = model
        self._boxes: List[tuple] = []

        outer = QVBoxLayout(self)
        intro = QLabel(
            f"<b>{mode.capitalize()} growth needs an organic carbon source.</b><br>"
            "Tick the substrate(s) to feed and set the uptake rate. Anything already "
            "open in the medium is ticked.")
        intro.setWordWrap(True)
        intro.setTextFormat(Qt.RichText)
        outer.addWidget(intro)

        if notes:
            applied = QLabel("<br>".join(f"• {n}" for n in notes))
            applied.setWordWrap(True)
            applied.setTextFormat(Qt.RichText)
            applied.setStyleSheet(f"color:{style.TEXT_MUTED};")
            outer.addWidget(applied)

        rate_row = QHBoxLayout()
        rate_row.addWidget(QLabel("Uptake rate (mmol gDW⁻¹ h⁻¹):"))
        self.rate = QDoubleSpinBox()
        self.rate.setDecimals(3)
        self.rate.setRange(0.001, 1000.0)
        self.rate.setValue(DEFAULT_UPTAKE)
        self.rate.setMaximumWidth(120)
        self.rate.setToolTip(
            "Applied to every ticked source.")
        rate_row.addWidget(self.rate)
        rate_row.addStretch(1)
        outer.addLayout(rate_row)

        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter…")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._filter)
        outer.addWidget(self.search)

        listing = QWidget()
        self._list_layout = QVBoxLayout(listing)
        self._list_layout.setContentsMargins(4, 4, 4, 4)
        for entry in self._candidates():
            box = QCheckBox(entry["label"])
            box.setChecked(entry["active"])
            box.setToolTip(f"{entry['id']} — {entry['carbon']} carbon atom(s)")
            self._list_layout.addWidget(box)
            self._boxes.append((box, entry))
        self._list_layout.addStretch(1)

        area = QScrollArea()
        area.setWidget(listing)
        area.setWidgetResizable(True)
        outer.addWidget(area, 1)

        self.count = QLabel("")
        self.count.setStyleSheet(f"color:{style.TEXT_MUTED};")
        outer.addWidget(self.count)
        for box, _ in self._boxes:
            box.toggled.connect(self._update_count)
        self._update_count()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Apply")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)
        clamp_to_screen(self)

    # -- data ------------------------------------------------------------------------
    def _candidates(self) -> List[dict]:
        """Organic carbon exchanges, the usual substrates first."""
        from ...core import physiology
        from ...core.network_graph import clean_label
        from ...core.physiology import _INORGANIC_C, _base as _met_base

        rows = []
        for row in physiology.substrate_exchanges(self._model):
            if not row.get("carbon"):
                continue
            # CO2 and bicarbonate carry carbon but are not *organic* carbon. Offering
            # them here would suggest feeding CO2 to a heterotroph, which is exactly the
            # confusion these modes exist to avoid — and it is the same test
            # `available_growth_modes` uses to decide the mode is possible at all.
            metabolite = next(iter(
                self._model.reactions.get_by_id(row["id"]).metabolites), None)
            if metabolite is not None and \
                    _met_base(metabolite.id).lower() in _INORGANIC_C:
                continue
            stem = row["id"]
            for prefix in ("EX_", "DM_", "SK_"):
                if stem.startswith(prefix):
                    stem = stem[len(prefix):]
                    break
            stem = stem[:-2] if stem.endswith("_e") else stem
            rows.append({
                "id": row["id"],
                "carbon": row["carbon"],
                "active": bool(row.get("active")),
                "common": stem in _COMMON,
                "label": clean_label(
                    f"{row['name']}  —  {row['id']}  ({row['carbon']} C)"),
            })
        rows.sort(key=lambda r: (not r["common"], not r["active"], r["label"].lower()))
        return rows

    def _filter(self, text: str) -> None:
        needle = text.strip().lower()
        for box, entry in self._boxes:
            box.setVisible(not needle or needle in entry["label"].lower())

    def _update_count(self) -> None:
        n = sum(1 for box, _ in self._boxes if box.isChecked())
        self.count.setText(f"{n} source(s) selected of {len(self._boxes)} available.")

    # -- result ----------------------------------------------------------------------
    def selection(self) -> Dict[str, float]:
        """``{exchange id: uptake rate}`` for every ticked source."""
        rate = self.rate.value()
        return {entry["id"]: rate for box, entry in self._boxes if box.isChecked()}

    @staticmethod
    def choose(parent, model, mode: str,
               notes: Sequence[str] = ()) -> Optional[Dict[str, float]]:
        dialog = CarbonSourceDialog(parent, model, mode, notes)
        if dialog.exec():
            return dialog.selection()
        return None
