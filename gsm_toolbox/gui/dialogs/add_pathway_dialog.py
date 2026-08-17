"""Confirmation + settings popup shown before a heterologous pathway is applied.

Summarises exactly what will change — the new reactions and new metabolites (whose
ids and names the user can edit right here), plus any compartment collapses / dropped
or unbalanced steps — and lets the user pick a category and choose how the product
leaves (or stays in) the cell: exported/secreted, volatile gas, accumulated
intracellularly, or no route.
"""

from __future__ import annotations

from typing import List

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QRadioButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ..widgets.dialog_util import clamp_to_screen

# (label, route-key, explanation) — offered in the Product export box.
_ROUTES = [
    ("Exported to the medium (secreted or diffuses out)", "secrete",
     "Adds a membrane transport step and an exchange so the product accumulates in the "
     "medium. The standard choice — needed for production envelopes, yields and strain "
     "design. Active secretion and passive diffusion are represented the same way."),
    ("Volatile — leaves the cell as a gas", "gas",
     "As export, but the exchange represents escape to the gas phase."),
    ("Accumulated inside the cell", "intracellular",
     "Adds an intracellular sink so the product can be made at steady state without "
     "leaving — for compounds stored inside (glycogen, PHB granules, pigments, storage "
     "lipids)."),
    ("No export route (advanced)", "none",
     "Adds nothing. Only valid if another reaction you add consumes the product; "
     "otherwise FBA cannot carry any flux to it."),
]


class AddPathwayDialog(QDialog):
    def __init__(self, parent, *, target_name: str, reactions: List[dict],
                 new_metabolites: List[dict], categories: List[str],
                 collapsed: List[str], dropped: List[str], unbalanced: List[tuple],
                 can_secrete: bool, product_exchange=None,
                 unverified: List[tuple] = None):
        super().__init__(parent)
        self.setWindowTitle("Add pathway to model")
        self.resize(680, 660)
        self._rxn_edits = []      # (orig_id, id_edit, name_edit)
        self._met_edits = []      # (orig_id, id_edit, name_edit)

        body = QWidget()
        v = QVBoxLayout(body)
        v.addWidget(QLabel(f"<b>Applying the pathway to {target_name}</b> will make the "
                           "following changes. You can edit any new id or name before adding."))

        # --- new reactions (editable id + name) ---
        rbox = QGroupBox(f"New reactions ({len(reactions)}) — edit id / name if you like")
        rg = QGridLayout(rbox)
        rg.addWidget(QLabel("<i>id</i>"), 0, 0)
        rg.addWidget(QLabel("<i>name</i>"), 0, 1)
        rg.addWidget(QLabel("<i>equation</i>"), 0, 2)
        for i, r in enumerate(reactions, start=1):
            oid = r["id"]
            id_edit = QLineEdit(oid); id_edit.setMaximumWidth(150)
            name_edit = QLineEdit(r.get("name", "") or "")
            eq = QLabel(r.get("equation", "")); eq.setWordWrap(True)
            eq.setStyleSheet("color:#5f6368;"); eq.setTextInteractionFlags(Qt.TextSelectableByMouse)
            rg.addWidget(id_edit, i, 0); rg.addWidget(name_edit, i, 1); rg.addWidget(eq, i, 2)
            self._rxn_edits.append((oid, id_edit, name_edit))
        rg.setColumnStretch(2, 1)
        v.addWidget(rbox)

        # --- new metabolites (editable id + name) ---
        mbox = QGroupBox(f"New metabolites ({len(new_metabolites)}) — edit id / name if you like")
        mg = QGridLayout(mbox)
        if new_metabolites:
            mg.addWidget(QLabel("<i>id</i>"), 0, 0)
            mg.addWidget(QLabel("<i>name</i>"), 0, 1)
            for i, m in enumerate(new_metabolites, start=1):
                oid = m["id"]
                id_edit = QLineEdit(oid); id_edit.setMaximumWidth(180)
                name_edit = QLineEdit(m.get("name", "") or "")
                mg.addWidget(id_edit, i, 0); mg.addWidget(name_edit, i, 1)
                self._met_edits.append((oid, id_edit, name_edit))
            mg.setColumnStretch(1, 1)
        else:
            mg.addWidget(QLabel("None — every metabolite already exists in the model."), 0, 0)
        v.addWidget(mbox)

        # --- other changes / warnings ---
        notes = []
        if collapsed:
            notes.append("Compartments collapsed into the host's main compartment: "
                         + ", ".join(collapsed))
        if dropped:
            notes.append(f"{len(dropped)} reaction(s) dropped (pure transport after collapse).")
        if unbalanced:
            notes.append("Some added reactions appear mass/charge-unbalanced — a suggestion to "
                         "check them: " + ", ".join(f"{rid} ({res})" for rid, res in unbalanced[:6]))
        # Separate note, and no ⚠: an unverifiable reaction is a missing formula in the
        # database, not a fault in the reaction. Merging the two turned every route into
        # a warning, which is how a real imbalance gets scrolled past.
        info_notes = []
        for rid, why in (unverified or [])[:6]:
            info_notes.append(f"{rid} cannot be balance-checked — {why}.")
        if notes or info_notes:
            nbox = QGroupBox("Other changes to review")
            nv = QVBoxLayout(nbox)
            for n in notes:
                lb = QLabel("⚠ " + n)
                lb.setWordWrap(True)
                nv.addWidget(lb)
            for n in info_notes:
                lb = QLabel("ℹ " + n)
                lb.setWordWrap(True)
                lb.setStyleSheet("color:#5f6368;")
                nv.addWidget(lb)
            v.addWidget(nbox)

        # --- category ---
        cbox = QGroupBox("Group the added reactions in a category")
        cv = QVBoxLayout(cbox)
        self.category_combo = QComboBox()
        self.category_combo.setEditable(True)
        default_cat = f"{target_name} production"
        self.category_combo.addItem(default_cat)
        for c in categories:
            if c != default_cat:
                self.category_combo.addItem(c)
        cv.addWidget(QLabel("Category name (pick an existing one or type a new name):"))
        cv.addWidget(self.category_combo)
        v.addWidget(cbox)

        # --- product export route ---
        ebox = QGroupBox("Product export")
        ev = QVBoxLayout(ebox)
        ev.addWidget(QLabel("How does the product leave (or stay in) the cell? "
                            "By default it is exported so you can measure a titre."))
        self._route_group = QButtonGroup(self)
        for label, key, desc in _ROUTES:
            rb = QRadioButton(label)
            rb.setProperty("route", key)
            self._route_group.addButton(rb)
            ev.addWidget(rb)
            d = QLabel(desc)
            d.setWordWrap(True)
            d.setStyleSheet("color:#5f6368; margin-left:22px; margin-bottom:4px;")
            ev.addWidget(d)
            if key == "secrete":
                rb.setChecked(True)         # default: product leaves the cell
        if product_exchange and product_exchange.get("id") and product_exchange.get("created"):
            ev.addWidget(QLabel(f"<span style='color:#5f6368'>Default export reaction: "
                                f"<b>{product_exchange['id']}</b> (created automatically).</span>"))
        elif product_exchange and product_exchange.get("id"):
            ev.addWidget(QLabel(f"<span style='color:#5f6368'>An exchange already exists and "
                                f"will be reused: <b>{product_exchange['id']}</b>.</span>"))
        v.addWidget(ebox)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(body)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Add pathway")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        outer = QVBoxLayout(self)
        outer.addWidget(scroll, 1)
        outer.addWidget(buttons)
        clamp_to_screen(self)

    def _route(self) -> str:
        btn = self._route_group.checkedButton()
        return btn.property("route") if btn else "secrete"

    @staticmethod
    def _collect(edits):
        renames, names = {}, {}
        for orig_id, id_edit, name_edit in edits:
            new_id = id_edit.text().strip()
            if new_id and new_id != orig_id:
                renames[orig_id] = new_id
            new_name = name_edit.text().strip()
            names[orig_id] = new_name   # keyed by original id; applied post-add
        return renames, names

    def values(self) -> dict:
        route = self._route()
        rxn_renames, rxn_names = self._collect(self._rxn_edits)
        met_renames, met_names = self._collect(self._met_edits)
        return {
            "category": self.category_combo.currentText().strip(),
            "route": route,
            "gas": route == "gas",             # kept for backward compatibility
            "rxn_renames": rxn_renames, "rxn_names": rxn_names,
            "met_renames": met_renames, "met_names": met_names,
        }
