"""Physiology and concentration settings for the MDF calculation.

The Max-min Driving Force is not a property of a pathway alone — it is a property of a
pathway *under assumed conditions*. pH, ionic strength and above all the allowed
metabolite concentration range change the answer, and widening the range always makes a
route look better. Burying those assumptions in a footnote invites a user to read an MDF
as if it were measured. This dialog puts them in front of the calculation: the settings
are shown first, and only a deliberate second click runs the analysis.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)

from .. import style
from ...core import thermodynamics as thermo

#: Defaults, matching the values the analysis used before it was configurable.
_DEFAULTS = {
    "ph": 7.5,
    "ionic_strength_M": 0.25,
    "pmg": 3.0,
    "temperature": thermo.T_DEFAULT,
    "conc_lo_mM": thermo.DEFAULT_CONC[0] * 1e3,      # 1e-6 M  → 0.001 mM
    "conc_hi_mM": thermo.DEFAULT_CONC[1] * 1e3,      # 1e-2 M  → 10 mM
}


class MDFSettingsDialog(QDialog):
    """Collects MDF assumptions. ``exec()`` returning Accepted means *run it now*."""

    def __init__(self, parent, *, metabolites: Optional[List[Tuple[str, str]]] = None,
                 initial: Optional[dict] = None):
        """``initial`` is a previous :meth:`settings` result, so a re-run starts where the
        last one left off rather than silently reverting to the defaults."""
        super().__init__(parent)
        self.setWindowTitle("Thermodynamic analysis — settings")
        self.resize(760, 620)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)
        values = dict(_DEFAULTS)
        values.update(_from_settings(initial))

        v = QVBoxLayout(self)
        v.addWidget(_wrap(
            "<b>These assumptions decide the answer.</b> The MDF is the driving force of "
            "the least-favourable step once concentrations have been optimised within the "
            "range you allow here — so a wider range always yields a better-looking "
            "number. Check them, then run the analysis."))

        # ---- physiology ----------------------------------------------------------
        phys = QGroupBox("Assumed physiology")
        form = QFormLayout(phys)
        self.ph = _spin(0.0, 14.0, 0.1, values["ph"], 2)
        self.ph.setToolTip("Cytosolic pH. eQuilibrator transforms ΔrG° to ΔrG′° at this pH.")
        form.addRow("pH", self.ph)
        self.ionic = _spin(0.0, 1.0, 0.01, values["ionic_strength_M"], 3)
        self.ionic.setToolTip("Ionic strength in mol/L — affects activity coefficients.")
        form.addRow("Ionic strength (M)", self.ionic)
        self.pmg = _spin(0.0, 10.0, 0.1, values["pmg"], 2)
        self.pmg.setToolTip("−log₁₀[Mg²⁺]. Matters most for phosphate-transfer reactions.")
        form.addRow("pMg", self.pmg)
        self.temp = _spin(273.15, 373.15, 0.5, values["temperature"], 2)
        self.temp.setToolTip("Temperature in kelvin (298.15 K = 25 °C).")
        form.addRow("Temperature (K)", self.temp)
        v.addWidget(phys)

        # ---- concentration range -------------------------------------------------
        conc = QGroupBox("Allowed metabolite concentrations")
        cform = QFormLayout(conc)
        self.conc_lo = _spin(1e-6, 1e4, 0.001, values["conc_lo_mM"], 6)
        self.conc_hi = _spin(1e-6, 1e4, 0.1, values["conc_hi_mM"], 6)
        self.conc_lo.setToolTip("Lower bound applied to every metabolite without an "
                                "explicit override below.")
        self.conc_hi.setToolTip("Upper bound applied to every metabolite without an "
                                "explicit override below.")
        cform.addRow("Minimum (mM)", self.conc_lo)
        cform.addRow("Maximum (mM)", self.conc_hi)
        cform.addRow(_wrap(
            "The default 0.001 – 10 mM spans the physiological range measured across "
            "central metabolism. Narrow it to test a route under realistic conditions; "
            "widen it only if you can justify the concentrations.", muted=True))
        v.addWidget(conc)

        # ---- per-metabolite overrides -------------------------------------------
        self._mets = list(metabolites or [])
        if self._mets:
            self.override_box = QCheckBox(
                f"Set concentrations for individual metabolites ({len(self._mets)} in "
                "this route)")
            self.override_box.setToolTip(
                "Pin a compound whose concentration you know — a fixed substrate feed, or "
                "a toxic intermediate that cannot be allowed to accumulate.")
            self.override_box.toggled.connect(self._toggle_overrides)
            v.addWidget(self.override_box)

            self.table = QTableWidget(len(self._mets), 3)
            self.table.setHorizontalHeaderLabels(
                ["Metabolite", "Minimum (mM)", "Maximum (mM)"])
            self.table.verticalHeader().setVisible(False)
            for r, (mid, name) in enumerate(self._mets):
                item = QTableWidgetItem(name or mid)
                item.setToolTip(mid)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(r, 0, item)
                for c in (1, 2):
                    blank = QTableWidgetItem("")
                    blank.setToolTip("Leave empty to use the default range above.")
                    self.table.setItem(r, c, blank)
            hh = self.table.horizontalHeader()
            hh.setSectionResizeMode(0, QHeaderView.Stretch)
            hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
            hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
            self.table.setVisible(False)
            v.addWidget(self.table, 1)
        else:
            self.override_box = None
            self.table = None
            v.addStretch(1)

        # ---- actions -------------------------------------------------------------
        row = QHBoxLayout()
        reset = QPushButton("Restore defaults")
        reset.setToolTip("Return every setting to the standard physiological assumption.")
        reset.clicked.connect(self._restore_defaults)
        row.addWidget(reset)
        row.addStretch(1)
        bb = QDialogButtonBox()
        self.run_btn = bb.addButton("Run analysis", QDialogButtonBox.AcceptRole)
        self.run_btn.setObjectName("primary")
        bb.addButton(QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        row.addWidget(bb)
        v.addLayout(row)

    # -- behaviour ---------------------------------------------------------------
    def _toggle_overrides(self, on: bool) -> None:
        if self.table is not None:
            self.table.setVisible(on)

    def _restore_defaults(self) -> None:
        self.ph.setValue(_DEFAULTS["ph"])
        self.ionic.setValue(_DEFAULTS["ionic_strength_M"])
        self.pmg.setValue(_DEFAULTS["pmg"])
        self.temp.setValue(_DEFAULTS["temperature"])
        self.conc_lo.setValue(_DEFAULTS["conc_lo_mM"])
        self.conc_hi.setValue(_DEFAULTS["conc_hi_mM"])
        if self.table is not None:
            for r in range(self.table.rowCount()):
                for c in (1, 2):
                    self.table.item(r, c).setText("")

    def accept(self) -> None:
        if self.conc_hi.value() <= self.conc_lo.value():
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Concentration range", "The maximum concentration must be greater "
                "than the minimum.")
            return
        super().accept()

    # -- results -----------------------------------------------------------------
    def settings(self) -> dict:
        """Keyword arguments for :func:`thermodynamics.analyse_pathway_mdf`.

        Concentrations are entered in mM because that is how they are measured and
        reported; the solver works in M, so they are converted here — in one place.
        """
        out = {
            "ph": float(self.ph.value()),
            "ionic_strength_M": float(self.ionic.value()),
            "pmg": float(self.pmg.value()),
            "temperature": float(self.temp.value()),
            "default_conc": (self.conc_lo.value() * 1e-3, self.conc_hi.value() * 1e-3),
        }
        bounds = self._overrides()
        if bounds:
            out["conc_bounds"] = bounds
        return out

    def _overrides(self) -> Dict[str, Tuple[float, float]]:
        if self.table is None or not (self.override_box and self.override_box.isChecked()):
            return {}
        lo_default = self.conc_lo.value() * 1e-3
        hi_default = self.conc_hi.value() * 1e-3
        out: Dict[str, Tuple[float, float]] = {}
        for r, (mid, _name) in enumerate(self._mets):
            lo = _as_float(self.table.item(r, 1))
            hi = _as_float(self.table.item(r, 2))
            if lo is None and hi is None:
                continue        # untouched row: the default range already covers it
            lo_m = lo * 1e-3 if lo is not None else lo_default
            hi_m = hi * 1e-3 if hi is not None else hi_default
            if hi_m < lo_m:
                lo_m, hi_m = hi_m, lo_m
            out[mid] = (lo_m, hi_m)
        return out

    def summary(self) -> str:
        """One line describing the assumptions, for the result window."""
        s = self.settings()
        lo, hi = s["default_conc"]
        extra = ""
        n = len(s.get("conc_bounds") or {})
        if n:
            extra = f", {n} metabolite(s) with individual bounds"
        return (f"pH {s['ph']:g}, ionic strength {s['ionic_strength_M']:g} M, "
                f"pMg {s['pmg']:g}, {s['temperature']:g} K, concentrations "
                f"{lo * 1e3:g}–{hi * 1e3:g} mM{extra}.")


# -- small helpers ---------------------------------------------------------------------
def _from_settings(settings: Optional[dict]) -> dict:
    """Turn a :meth:`MDFSettingsDialog.settings` result back into widget values.

    The two shapes differ deliberately: the analysis takes a concentration range in M as
    one tuple, while the form edits a minimum and a maximum in mM.
    """
    if not settings:
        return {}
    out = {k: settings[k] for k in ("ph", "ionic_strength_M", "pmg", "temperature")
           if k in settings}
    conc = settings.get("default_conc")
    if conc and len(conc) == 2:
        out["conc_lo_mM"] = float(conc[0]) * 1e3
        out["conc_hi_mM"] = float(conc[1]) * 1e3
    return out


def _spin(lo: float, hi: float, step: float, value: float, decimals: int) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setDecimals(decimals)
    s.setRange(lo, hi)
    s.setSingleStep(step)
    s.setValue(value)
    return s


def _wrap(html: str, *, muted: bool = False) -> QLabel:
    lbl = QLabel(html)
    lbl.setWordWrap(True)
    lbl.setTextFormat(Qt.RichText)
    if muted:
        lbl.setStyleSheet(f"color:{style.TEXT_MUTED};")
    return lbl


def _as_float(item) -> Optional[float]:
    text = (item.text() if item is not None else "").strip()
    if not text:
        return None
    try:
        val = float(text.replace(",", "."))
    except ValueError:
        return None
    return val if val > 0 else None
