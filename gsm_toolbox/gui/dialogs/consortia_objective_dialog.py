"""Objective settings for a consortium (community) model.

A plain community FBA maximises the *summed* biomass of all members, which lets the
faster grower take all the resources and starve the others. This dialog exposes the
two controls an expert uses to model a consortium realistically:

* **Dominance weight** per member — how strongly that organism's growth counts in
  the community objective (relative abundance). Equal weights = balanced; a higher
  weight makes that species dominate.
* **Minimum growth per member** — each organism must grow at least this fraction of
  the maximum it could reach in the community, so no member is starved (a
  SteadyCom-style floor).
"""

from __future__ import annotations

from typing import Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QVBoxLayout,
)


class ConsortiaObjectiveDialog(QDialog):
    def __init__(self, parent, members: List[str], *,
                 weights: Optional[Dict[str, float]] = None,
                 min_growth_fraction: float = 0.0):
        super().__init__(parent)
        self.setWindowTitle("Consortia objective settings")
        self.resize(460, 420)
        weights = weights or {}

        v = QVBoxLayout(self)
        intro = QLabel(
            "A community FBA maximises the members' combined growth. Set each "
            "organism's <b>dominance</b> (relative abundance in the community "
            "objective) and a <b>minimum growth</b> so no member is starved.")
        intro.setWordWrap(True)
        v.addWidget(intro)

        wbox = QGroupBox("Dominance weight per member")
        wform = QFormLayout(wbox)
        self._weight_spins: Dict[str, QDoubleSpinBox] = {}
        for name in members:
            sp = QDoubleSpinBox()
            sp.setRange(0.0, 100.0)
            sp.setDecimals(2)
            sp.setSingleStep(0.5)
            sp.setValue(float(weights.get(name, 1.0)))
            sp.setToolTip("Relative weight of this organism's growth in the community "
                          "objective. Equal values = balanced; higher = more dominant.")
            self._weight_spins[name] = sp
            wform.addRow(f"{name}:", sp)
        v.addWidget(wbox)

        mbox = QGroupBox("No-starvation floor")
        mform = QFormLayout(mbox)
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(0.0, 0.95)
        self.min_spin.setDecimals(2)
        self.min_spin.setSingleStep(0.05)
        self.min_spin.setValue(float(min_growth_fraction or 0.0))
        self.min_spin.setToolTip("Each member must grow at least this fraction of the "
                                 "maximum it could reach in the community (0 = no floor).")
        mform.addRow("Minimum growth per member\n(fraction of its own max):", self.min_spin)
        note = QLabel("Tip: if the community becomes infeasible, the floor is too high "
                      "for the shared nutrients — lower it.")
        note.setWordWrap(True)
        note.setStyleSheet("color:#6B7280; font-size:11px;")
        mform.addRow(note)
        v.addWidget(mbox)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Apply")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)

    def values(self) -> dict:
        return {
            "weights": {name: sp.value() for name, sp in self._weight_spins.items()},
            "min_growth_fraction": self.min_spin.value(),
        }
