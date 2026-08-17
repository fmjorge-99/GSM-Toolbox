"""Interactive reaction builder.

The user assembles a reaction from dropdown rows (coefficient + searchable
metabolite picker on the reactant and product sides), chooses reversibility,
bounds, gene rule and subsystem, and sees a live preview with the 2-D structures
of the metabolites. Metabolites can be searched by name in the model, a loaded
reaction database, or KEGG; a KEGG hit keeps the *compound name* and suggests an
editable metabolite id, recording the KEGG id as an annotation. A reaction can
also be loaded automatically from an EC number.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import cobra
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QComboBox,
    QCompleter,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ...core import namespace
from ..widgets.structure_fetcher import StructureFetcher


def _display(met_id: str, name: str) -> str:
    from ...core.network_graph import short_metabolite_name
    base = met_id.rsplit("_", 1)[0]
    name = short_metabolite_name(met_id, name)
    if name and name not in (met_id, base):
        return f"{name}  ({met_id})"
    return met_id


def _slug_id(name: str, compartment: str = "c") -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", (name or "").lower()).strip("_")[:24]
    return f"{slug or 'met'}_{compartment}"


def _parse_formula(formula: str) -> Dict[str, float]:
    """Element -> count from a chemical formula (e.g. 'C6H12O6' -> {C:6,H:12,O:6})."""
    counts: Dict[str, float] = {}
    for elem, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        if not elem:
            continue
        counts[elem] = counts.get(elem, 0.0) + (float(num) if num else 1.0)
    return counts


class _MetRow(QWidget):
    """One term: coefficient spin box + searchable metabolite combo + remove."""

    def __init__(self, palette: List[Tuple[str, str]], on_remove, on_change):
        super().__init__()
        self.coeff = QDoubleSpinBox()
        self.coeff.setRange(0.0001, 100000.0)
        self.coeff.setDecimals(4)
        self.coeff.setValue(1.0)
        self.coeff.setMaximumWidth(90)
        self.coeff.valueChanged.connect(lambda _: on_change())

        self.combo = QComboBox()
        self.combo.setEditable(True)
        self.combo.setInsertPolicy(QComboBox.NoInsert)
        self.combo.setMinimumWidth(260)
        for display, mid in palette:
            self.combo.addItem(display, mid)
        self.combo.setCurrentIndex(-1)
        self.combo.setCurrentText("")
        completer = self.combo.completer()
        if completer is not None:
            completer.setCompletionMode(QCompleter.PopupCompletion)
            completer.setFilterMode(Qt.MatchContains)
        self.combo.currentTextChanged.connect(lambda _: on_change())

        remove = QToolButton()
        remove.setText("✕")
        remove.setToolTip("Remove this metabolite")
        remove.clicked.connect(lambda: on_remove(self))

        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.coeff)
        lay.addWidget(self.combo, 1)
        lay.addWidget(remove)

    def metabolite_id(self) -> str:
        idx = self.combo.currentIndex()
        text = self.combo.currentText().strip()
        if idx >= 0 and self.combo.itemText(idx) == text:
            data = self.combo.itemData(idx)
            if data:
                return str(data)
        if text.endswith(")") and "(" in text:
            return text[text.rfind("(") + 1:-1].strip()
        return text

    def set_metabolite(self, display: str, mid: str) -> None:
        if self.combo.findData(mid) < 0:
            self.combo.addItem(display, mid)
        self.combo.setCurrentIndex(self.combo.findData(mid))

    def add_option(self, display: str, mid: str) -> None:
        if self.combo.findData(mid) < 0:
            self.combo.addItem(display, mid)


class ReactionBuilderDialog(QDialog):
    """Build a new reaction from metabolite dropdowns and stoichiometry."""

    def __init__(self, model: cobra.Model, parent=None,
                 database: Optional[cobra.Model] = None):
        super().__init__(parent)
        self.setWindowTitle("Add Reaction — builder")
        self.setMinimumSize(760, 660)
        self._model = model
        self._database = database
        self._palette = self._build_palette(model)
        self._reactant_rows: List[_MetRow] = []
        self._product_rows: List[_MetRow] = []
        self._ec = ""
        # metadata for metabolites the builder introduces (name + annotation),
        # applied to newly-created metabolites when the reaction is added.
        self._met_info: Dict[str, dict] = {}
        # cross-ref index of the host: KEGG id / inchikey -> existing metabolite id
        self._host_xref = self._build_host_xref(model)
        self._struct_cache: Dict[str, QPixmap] = {}
        self._fetchers: List[_StructureFetcher] = []

        self.id_edit = QLineEdit()
        self.id_edit.setPlaceholderText("Reaction id")
        self.name_edit = QLineEdit()
        self.id_edit.textChanged.connect(self._update_preview)

        meta = QFormLayout()
        meta.addRow("Reaction ID:", self.id_edit)
        meta.addRow("Name:", self.name_edit)

        # action buttons
        search_btn = QPushButton("Search a metabolite (database / KEGG)…")
        search_btn.setToolTip("Find a metabolite by name in a loaded reaction database or "
                              "from KEGG, and add it to the pickers below.")
        search_btn.clicked.connect(self._search_metabolite)
        ec_btn = QPushButton("Load reaction by EC number…")
        ec_btn.setToolTip("Paste an EC number; the reaction it catalyses is fetched "
                          "and loaded into the builder.")
        ec_btn.clicked.connect(self._load_by_ec)
        actions = QHBoxLayout()
        actions.addWidget(search_btn)
        actions.addWidget(ec_btn)

        # reactant / product columns
        self._reactant_box = self._side_box("Reactants (consumed)", self._reactant_rows,
                                            self._add_reactant_row)
        self._product_box = self._side_box("Products (produced)", self._product_rows,
                                           self._add_product_row)
        arrow = QLabel("→")
        arrow.setStyleSheet("font-size:28px; color:#1597B8;")
        arrow.setAlignment(Qt.AlignCenter)
        cols = QHBoxLayout()
        cols.addWidget(self._reactant_box, 1)
        cols.addWidget(arrow)
        cols.addWidget(self._product_box, 1)

        # direction + bounds
        self.rev_reversible = QRadioButton("Reversible  (↔)")
        self.rev_irreversible = QRadioButton("Irreversible  (→)")
        self.rev_irreversible.setChecked(True)
        self.rev_reversible.toggled.connect(self._on_direction_changed)
        self.lower = QDoubleSpinBox()
        self.lower.setRange(-1e6, 1e6)
        self.lower.setDecimals(2)
        self.lower.setValue(0.0)
        self.upper = QDoubleSpinBox()
        self.upper.setRange(-1e6, 1e6)
        self.upper.setDecimals(2)
        self.upper.setValue(1000.0)
        dir_row = QHBoxLayout()
        dir_row.addWidget(self.rev_irreversible)
        dir_row.addWidget(self.rev_reversible)
        dir_row.addStretch(1)
        dir_row.addWidget(QLabel("Lower:"))
        dir_row.addWidget(self.lower)
        dir_row.addWidget(QLabel("Upper:"))
        dir_row.addWidget(self.upper)

        # gene rule + subsystem
        self.gpr_edit = QLineEdit()
        self.gpr_edit.setPlaceholderText("e.g. b0001 and b0002")
        self.subsystem_edit = QLineEdit()
        meta2 = QFormLayout()
        meta2.addRow("Gene rule (GPR):", self.gpr_edit)
        meta2.addRow("Subsystem:", self.subsystem_edit)
        self.ec_label = QLabel("EC number: —")
        self.ec_label.setStyleSheet("color:#5f6368;")

        # live text preview + on-demand structure view
        self.preview = QLabel()
        self.preview.setWordWrap(True)
        self.preview.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.preview.setStyleSheet(
            "background:#F2F5FB; border:1px solid #D6E0F0; border-radius:6px; padding:8px;")
        self.show_rxn_btn = QPushButton("Show reaction")
        self.show_rxn_btn.setToolTip("Fetch the metabolites' 2-D structures and draw the "
                                     "reaction (structures are cached for offline reuse).")
        self.show_rxn_btn.clicked.connect(self._refresh_structures)
        preview_row = QHBoxLayout()
        preview_row.addWidget(QLabel("Preview:"))
        preview_row.addStretch(1)
        preview_row.addWidget(self.show_rxn_btn)

        # Live mass/charge balance indicator (Issue 12) + one-click fixer.
        self.balance_label = QLabel("Balance: —")
        self.balance_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.balance_btn = QPushButton("Balance with H⁺/H₂O")
        self.balance_btn.setToolTip("Add or adjust water and protons to close a hydrogen/"
                                    "oxygen/charge imbalance.")
        self.balance_btn.clicked.connect(self._balance_with_water_proton)
        balance_row = QHBoxLayout()
        balance_row.addWidget(self.balance_label, 1)
        balance_row.addWidget(self.balance_btn)

        self._struct_host = QWidget()
        self._struct_layout = QHBoxLayout(self._struct_host)
        self._struct_layout.setContentsMargins(4, 4, 4, 4)
        self._struct_layout.addWidget(QLabel("Click “Show reaction” to render the structures."))
        struct_scroll = QScrollArea()
        struct_scroll.setWidgetResizable(True)
        struct_scroll.setFixedHeight(180)
        struct_scroll.setWidget(self._struct_host)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)

        # All the form content goes in a scroll area; the OK/Cancel buttons stay
        # pinned below it so they're always reachable even on short screens (#B6).
        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.addLayout(meta)
        content_layout.addLayout(actions)
        content_layout.addLayout(cols, 1)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        content_layout.addWidget(line)
        content_layout.addLayout(dir_row)
        content_layout.addLayout(meta2)
        content_layout.addWidget(self.ec_label)
        content_layout.addLayout(preview_row)
        content_layout.addWidget(self.preview)
        content_layout.addLayout(balance_row)
        content_layout.addWidget(struct_scroll)

        content_scroll = QScrollArea()
        content_scroll.setWidgetResizable(True)
        content_scroll.setFrameShape(QScrollArea.NoFrame)
        content_scroll.setWidget(content)

        layout = QVBoxLayout(self)
        layout.addWidget(content_scroll, 1)
        layout.addWidget(buttons)

        # Cap the initial size to the screen so OK/Cancel are never off-screen.
        from ..widgets.dialog_util import clamp_to_screen
        self.resize(720, 720)
        clamp_to_screen(self)

        self._add_reactant_row()
        self._add_product_row()
        self._update_preview()

    # ----- palette / host xref ----------------------------------------
    def _build_palette(self, model: cobra.Model) -> List[Tuple[str, str]]:
        items = [(_display(m.id, m.name or ""), m.id) for m in model.metabolites]
        items.sort(key=lambda it: it[0].lower())
        return items

    def _build_host_xref(self, model: cobra.Model) -> Dict[str, str]:
        index: Dict[str, str] = {}
        for m in model.metabolites:
            for tok in namespace.metabolite_tokens(m):
                index.setdefault(tok, m.id)
        return index

    def _all_rows(self) -> List[_MetRow]:
        return self._reactant_rows + self._product_rows

    # ----- side columns ------------------------------------------------
    def _side_box(self, title: str, store: List[_MetRow], add_cb) -> QGroupBox:
        box = QGroupBox(title)
        v = QVBoxLayout(box)
        host = QWidget()
        rows_layout = QVBoxLayout(host)
        rows_layout.setContentsMargins(0, 0, 0, 0)
        rows_layout.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(host)
        scroll.setFrameShape(QScrollArea.NoFrame)
        v.addWidget(scroll, 1)
        add = QPushButton("+ Add metabolite")
        add.clicked.connect(add_cb)
        v.addWidget(add)
        if "Reactant" in title:
            self._reactant_layout = rows_layout
        else:
            self._product_layout = rows_layout
        return box

    def _add_reactant_row(self) -> _MetRow:
        return self._add_row(self._reactant_rows, self._reactant_layout)

    def _add_product_row(self) -> _MetRow:
        return self._add_row(self._product_rows, self._product_layout)

    def _add_row(self, store: List[_MetRow], layout) -> _MetRow:
        row = _MetRow(self._palette, self._remove_row, self._update_preview)
        store.append(row)
        layout.insertWidget(layout.count() - 1, row)  # before the stretch
        self._update_preview()
        return row

    def _remove_row(self, row: _MetRow) -> None:
        for store in (self._reactant_rows, self._product_rows):
            if row in store:
                store.remove(row)
        row.setParent(None)
        row.deleteLater()
        self._update_preview()

    # ----- direction ---------------------------------------------------
    def _on_direction_changed(self) -> None:
        if self.rev_reversible.isChecked():
            self.lower.setValue(-1000.0)
            self.upper.setValue(1000.0)
        else:
            self.lower.setValue(0.0)
            self.upper.setValue(1000.0)
        self._update_preview()

    # ----- register a metabolite (resolve to host, or create new id) ---
    def _register_metabolite(self, name: str, kegg_id: str = "", inchikey: str = "",
                             suggest_id: str = "", ask: bool = True) -> Tuple[str, str]:
        """Return (display, met_id). If the compound matches an existing model
        metabolite (via KEGG/InChIKey xref) that id is reused so the reaction
        connects; otherwise a name-based id is suggested (and made editable)."""
        # try to connect to an existing model metabolite
        for tok in ([f"kegg:{kegg_id}"] if kegg_id else []) + \
                   ([f"inchikey:{inchikey.upper()}"] if inchikey else []):
            if tok in self._host_xref:
                mid = self._host_xref[tok]
                m = self._model.metabolites.get_by_id(mid)
                display = _display(mid, m.name or name)
                return display, mid

        mid = suggest_id or _slug_id(name)
        if ask:
            mid, ok = QInputDialog.getText(
                self, "Metabolite id",
                f"Metabolite id for “{name}” (edit if needed):", text=mid)
            if not ok or not mid.strip():
                return "", ""
            mid = mid.strip()
        ann: Dict[str, list] = {}
        if kegg_id:
            ann["kegg.compound"] = [kegg_id]
        if inchikey:
            ann["inchi_key"] = [inchikey]
        self._met_info[mid] = {"name": name, "annotation": ann}
        display = _display(mid, name)
        if not any(p[1] == mid for p in self._palette):
            self._palette.append((display, mid))
        return display, mid

    # ----- metabolite search ------------------------------------------
    def _search_metabolite(self) -> None:
        term, ok = QInputDialog.getText(
            self, "Search metabolite",
            "Metabolite name or id (searches a loaded database, then KEGG):")
        if not ok or not term.strip():
            return
        term = term.strip()
        low = term.lower()
        # (label, kind, payload) — kind: 'db' -> met object; 'kegg' -> (cid, name)
        matches = []
        if self._database is not None:
            for m in self._database.metabolites:
                if low in m.id.lower() or low in (m.name or "").lower():
                    matches.append((f"{m.name or m.id}  ({m.id})", "db", m))
                if len(matches) >= 40:
                    break
        if not matches:
            try:
                from ...core import databases
                for cid, desc in databases.kegg_find_compound(term)[:25]:
                    name = desc.split(";")[0].strip()
                    matches.append((f"{name}   [KEGG {cid}]", "kegg", (cid, name)))
            except Exception as exc:  # noqa: BLE001
                QMessageBox.warning(self, "Online search failed", str(exc))
                return
        if not matches:
            QMessageBox.information(self, "No matches",
                                    f"No metabolite matching '{term}' was found.")
            return
        labels = [m[0] for m in matches]
        choice, ok = QInputDialog.getItem(self, "Select metabolite",
                                          "Add which metabolite?", labels, 0, False)
        if not ok:
            return
        _, kind, payload = matches[labels.index(choice)]
        if kind == "db":
            m = payload
            display = _display(m.id, m.name or "")
            self._met_info[m.id] = {"name": m.name or m.id,
                                    "annotation": dict(getattr(m, "annotation", {}) or {})}
            if not any(p[1] == m.id for p in self._palette):
                self._palette.append((display, m.id))
            mid = m.id
        else:
            cid, name = payload
            display, mid = self._register_metabolite(name, kegg_id=cid)
            if not mid:
                return
        for row in self._all_rows():
            row.add_option(display, mid)
        QMessageBox.information(self, "Metabolite added",
                                f"“{display}” is now selectable in the pickers.")

    # ----- load reaction from an EC number ----------------------------
    def _load_by_ec(self) -> None:
        ec, ok = QInputDialog.getText(
            self, "Load reaction by EC number",
            "EC number (e.g. 4.1.1.102 or EC:4.1.1.102):")
        if not ok or not ec.strip():
            return
        from ...core import databases
        try:
            rids = databases.kegg_reactions_for_ec(ec)
            reactions = databases.kegg_get_reactions(rids[:6]) if rids else []
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "EC lookup failed", str(exc))
            return
        if not reactions:
            QMessageBox.information(
                self, "No reaction found",
                f"KEGG has no reaction linked to EC {ec.strip()}.")
            return
        if len(reactions) == 1:
            kr = reactions[0]
        else:
            labels = [f"{r.rid}: {r.name or r.equation}" for r in reactions]
            choice, ok = QInputDialog.getItem(
                self, "Choose reaction",
                f"EC {ec.strip()} catalyses several reactions — pick one:", labels, 0, False)
            if not ok:
                return
            kr = reactions[labels.index(choice)]
        self._populate_from_kegg(kr, ec.strip())

    def _populate_from_kegg(self, kr, ec: str) -> None:
        from ...core import databases
        cids = list(kr.substrates) + list(kr.products)
        try:
            names = databases.kegg_compound_names(cids)
        except Exception:  # noqa: BLE001
            names = {}
        # clear existing rows
        for row in list(self._all_rows()):
            self._remove_row(row)
        if not self.id_edit.text().strip():
            self.id_edit.setText(kr.rid)
        if not self.name_edit.text().strip():
            self.name_edit.setText(kr.name or kr.rid)
        self.rev_reversible.setChecked(bool(kr.reversible))
        self._ec = ec
        self.ec_label.setText(f"EC number: {ec}  (from KEGG {kr.rid})")

        def _fill(side: Dict[str, float], add_row):
            for cid, coeff in side.items():
                name = names.get(cid, cid)
                display, mid = self._register_metabolite(name, kegg_id=cid, ask=False)
                if not mid:
                    continue
                row = add_row()
                row.coeff.setValue(float(coeff))
                row.set_metabolite(display, mid)

        _fill(kr.substrates, self._add_reactant_row)
        _fill(kr.products, self._add_product_row)
        if not self._reactant_rows:
            self._add_reactant_row()
        if not self._product_rows:
            self._add_product_row()
        self._update_preview()

    # ----- preview + structures ---------------------------------------
    def _terms(self, rows: List[_MetRow]) -> List[Tuple[float, str]]:
        return [(row.coeff.value(), row.metabolite_id()) for row in rows
                if row.metabolite_id()]

    def _equation(self) -> str:
        arrow = "<=>" if self.rev_reversible.isChecked() else "-->"
        left = " + ".join(f"{c:g} {m}" for c, m in self._terms(self._reactant_rows))
        right = " + ".join(f"{c:g} {m}" for c, m in self._terms(self._product_rows))
        return f"{left} {arrow} {right}".strip()

    def _met_name(self, mid: str) -> str:
        if mid in self._met_info:
            return self._met_info[mid]["name"]
        if self._model.metabolites.has_id(mid):
            return self._model.metabolites.get_by_id(mid).name or mid
        return mid

    def _update_preview(self) -> None:
        rid = self.id_edit.text().strip() or "(reaction id)"
        self.preview.setText(f"<b>{rid}:</b>  {self._equation()}")
        self._update_balance()

    # ----- mass/charge balance (Issue 12) ------------------------------
    def _met_formula_charge(self, mid: str):
        """Return ``(elements_dict, charge, known)`` for a metabolite id."""
        if self._model.metabolites.has_id(mid):
            m = self._model.metabolites.get_by_id(mid)
            try:
                els = dict(m.elements)
            except Exception:  # noqa: BLE001
                els = {}
            if not els and not (m.formula or ""):
                return {}, 0.0, False
            return els, float(m.charge or 0), True
        info = self._met_info.get(mid, {})
        formula = info.get("formula")
        if formula:
            return _parse_formula(formula), float(info.get("charge", 0) or 0), True
        return {}, 0.0, False

    def _compute_balance(self):
        """Net element/charge residual (products − reactants) and any metabolites
        whose formula is unknown (so the verdict can't be trusted)."""
        net: Dict[str, float] = {}
        unknown: List[str] = []
        for sign, rows in ((-1.0, self._reactant_rows), (1.0, self._product_rows)):
            for coeff, mid in self._terms(rows):
                els, charge, known = self._met_formula_charge(mid)
                if not known:
                    unknown.append(mid)
                    continue
                for e, n in els.items():
                    net[e] = net.get(e, 0.0) + sign * coeff * n
                net["charge"] = net.get("charge", 0.0) + sign * coeff * charge
        net = {k: v for k, v in net.items() if abs(v) > 1e-6}
        return net, unknown

    def _update_balance(self) -> None:
        if not (self._terms(self._reactant_rows) or self._terms(self._product_rows)):
            self.balance_label.setText("Balance: —")
            self.balance_label.setStyleSheet("color:#5f6368;")
            self.balance_btn.setEnabled(False)
            return
        net, unknown = self._compute_balance()
        if unknown:
            # "Cannot be checked", never "unbalanced": the unknown participant's atoms
            # were never counted, so any residual would be that participant's own mass.
            names = ", ".join(sorted(set(self._met_name(m) for m in unknown))[:4])
            self.balance_label.setText(
                f"Balance: cannot be checked — no formula for {names}")
            self.balance_label.setStyleSheet("color:#B06000;")
            self.balance_btn.setEnabled(False)
            return
        if not net:
            self.balance_label.setText("Balance: ✓ mass and charge balanced")
            self.balance_label.setStyleSheet("color:#1E8E3E; font-weight:600;")
            self.balance_btn.setEnabled(False)
            return
        residual = ", ".join(f"{k}{v:+.3g}" for k, v in sorted(net.items()))
        self.balance_label.setText(f"Balance: ✗ unbalanced — {residual}")
        self.balance_label.setStyleSheet("color:#D93025; font-weight:600;")
        self.balance_btn.setEnabled(self._water_proton_fix(net) is not None)

    @staticmethod
    def _water_proton_fix(net: Dict[str, float]):
        """Return ``(n_water, n_proton)`` that would close ``net`` using only H₂O/H⁺,
        or ``None`` if the residual can't be fixed that way. Positive counts go on
        the product side, negative on the reactant side."""
        if set(net) - {"H", "O", "charge"}:
            return None
        n_w = -net.get("O", 0.0)
        n_h = -net.get("charge", 0.0)
        if abs(net.get("H", 0.0) + 2 * n_w + n_h) > 1e-6:
            return None
        if abs(n_w) < 1e-9 and abs(n_h) < 1e-9:
            return None
        return n_w, n_h

    def _dominant_compartment(self) -> str:
        from collections import Counter
        cnt = Counter()
        for _c, mid in self._terms(self._reactant_rows) + self._terms(self._product_rows):
            if self._model.metabolites.has_id(mid):
                cnt[self._model.metabolites.get_by_id(mid).compartment or "c"] += 1
        return cnt.most_common(1)[0][0] if cnt else "c"

    def _balance_with_water_proton(self) -> None:
        net, unknown = self._compute_balance()
        if unknown:
            return
        fix = self._water_proton_fix(net)
        if fix is None:
            QMessageBox.information(
                self, "Cannot auto-balance",
                "This residual cannot be closed with water and protons alone "
                "(other elements, or an inconsistent H/O/charge residual). "
                "Adjust the stoichiometry manually.")
            return
        comp = self._dominant_compartment()
        h2o, h = f"h2o_{comp}", f"h_{comp}"
        if not (self._model.metabolites.has_id(h2o) and self._model.metabolites.has_id(h)):
            QMessageBox.information(
                self, "Cannot auto-balance",
                f"The model has no water ({h2o}) and/or proton ({h}) to balance with.")
            return
        n_w, n_h = fix
        self._apply_balancing_term(h2o, n_w)   # water fixes oxygen (and part of H)
        self._apply_balancing_term(h, n_h)     # proton fixes charge (and remaining H)
        self._update_preview()

    def _apply_balancing_term(self, mid: str, amount: float) -> None:
        """Add ``amount`` of ``mid`` to the product side (or |amount| to the reactant
        side if negative), merging into an existing row for that metabolite."""
        if abs(amount) < 1e-9:
            return
        rows = self._product_rows if amount > 0 else self._reactant_rows
        add_row = self._add_product_row if amount > 0 else self._add_reactant_row
        for row in rows:
            if row.metabolite_id() == mid:
                row.coeff.setValue(row.coeff.value() + abs(amount))
                return
        display = mid
        for disp, pid in self._palette:
            if pid == mid:
                display = disp
                break
        row = add_row()
        row.coeff.setValue(abs(amount))
        row.set_metabolite(display, mid)

    def _refresh_structures(self) -> None:
        while self._struct_layout.count():
            item = self._struct_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        arrow_shown = False
        reactants = self._terms(self._reactant_rows)
        products = self._terms(self._product_rows)
        for i, (coeff, mid) in enumerate(reactants):
            self._add_structure_slot(coeff, mid)
            if i < len(reactants) - 1:
                self._struct_layout.addWidget(self._sign("+"))
        self._struct_layout.addWidget(self._sign("→"))
        arrow_shown = True
        for i, (coeff, mid) in enumerate(products):
            self._add_structure_slot(coeff, mid)
            if i < len(products) - 1:
                self._struct_layout.addWidget(self._sign("+"))
        self._struct_layout.addStretch(1)
        if not reactants and not products and not arrow_shown:
            self._struct_layout.addWidget(QLabel("Structures appear here as you add metabolites."))

    def _sign(self, text: str) -> QLabel:
        lab = QLabel(text)
        lab.setStyleSheet("font-size:20px; color:#5f6368;")
        lab.setAlignment(Qt.AlignCenter)
        return lab

    def _add_structure_slot(self, coeff: float, mid: str) -> None:
        cell = QWidget()
        v = QVBoxLayout(cell)
        v.setContentsMargins(2, 2, 2, 2)
        img = QLabel()
        img.setFixedSize(120, 120)
        img.setAlignment(Qt.AlignCenter)
        img.setStyleSheet("border:1px solid #E0E3E8; border-radius:6px; background:white;")
        name = self._met_name(mid)
        if mid in self._struct_cache:
            img.setPixmap(self._struct_cache[mid])
        else:
            img.setText("…")
            self._start_structure_fetch(mid, name)
        coeff_txt = "" if abs(coeff - 1.0) < 1e-9 else f"{coeff:g} × "
        cap = QLabel(f"{coeff_txt}{name}")
        cap.setWordWrap(True)
        cap.setAlignment(Qt.AlignCenter)
        cap.setMaximumWidth(130)
        cap.setStyleSheet("font-size:10px;")
        v.addWidget(img)
        v.addWidget(cap)
        cell.setProperty("met_img", True)
        img.setObjectName(f"img::{mid}")
        self._struct_layout.addWidget(cell)

    def _start_structure_fetch(self, mid: str, name: str) -> None:
        # Prefer the model metabolite's full annotation (SMILES/InChI/KEGG/ChEBI) so
        # the correct structure is fetched, not a name-collision hit (#B11).
        if self._model.metabolites.has_id(mid):
            fetcher = StructureFetcher.for_metabolite(
                self._model.metabolites.get_by_id(mid), size=140)
        else:
            info = self._met_info.get(mid) or {}
            ann = info.get("annotation") or {}

            def _a(*keys):
                for k in keys:
                    v = ann.get(k)
                    if v:
                        return v[0] if isinstance(v, (list, tuple)) else str(v)
                return ""
            fetcher = StructureFetcher(
                mid, name=name, size=140, inchikey=_a("inchi_key", "inchikey"),
                smiles=_a("smiles"), inchi=_a("inchi"), kegg=_a("kegg.compound", "kegg"))
        fetcher.fetched.connect(self._on_structure)
        self._fetchers.append(fetcher)
        fetcher.start()

    def _on_structure(self, mid: str, data: bytes) -> None:
        if not data:
            return
        pix = QPixmap()
        if not pix.loadFromData(data):
            return
        pix = pix.scaled(116, 116, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self._struct_cache[mid] = pix
        for img in self._struct_host.findChildren(QLabel):
            if img.objectName() == f"img::{mid}":
                img.setText("")
                img.setPixmap(pix)

    # ----- result ------------------------------------------------------
    def _on_accept(self) -> None:
        if not self.id_edit.text().strip():
            QMessageBox.warning(self, "Missing id", "Please enter a reaction ID.")
            return
        if not self._terms(self._reactant_rows) and not self._terms(self._product_rows):
            QMessageBox.warning(self, "Empty reaction",
                                "Add at least one metabolite to the reaction.")
            return
        self.accept()

    def values(self) -> dict:
        return {
            "reaction_id": self.id_edit.text().strip(),
            "name": self.name_edit.text().strip(),
            "reaction_string": self._equation(),
            "gene_reaction_rule": self.gpr_edit.text().strip(),
            "subsystem": self.subsystem_edit.text().strip(),
            "lower_bound": self.lower.value(),
            "upper_bound": self.upper.value(),
            "ec_number": self._ec,
            # name/annotation for metabolites the builder introduced (for A4: so
            # new metabolites are created with proper names + cross-references)
            "metabolite_info": dict(self._met_info),
        }
