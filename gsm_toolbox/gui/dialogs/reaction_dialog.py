"""Add / edit a reaction via a plain-language form.

Returns a dict of field values; the caller applies the change through the core
``editing`` module inside a ``Project.apply_edit`` so it is undoable.
"""

from __future__ import annotations

from typing import Optional

import cobra
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from ...core import balancing


class ReactionDialog(QDialog):
    """Create a new reaction or edit an existing one."""

    def __init__(self, parent=None, reaction: Optional[cobra.Reaction] = None):
        super().__init__(parent)
        self._editing = reaction is not None
        self._reaction = reaction
        self._balance_changes = []      # staged H+/H2O additions (applied on OK)
        self.setWindowTitle("Edit Reaction" if self._editing else "Add Reaction")
        self.setMinimumWidth(480)

        self.id_edit = QLineEdit()
        self.name_edit = QLineEdit()
        self.equation_edit = QLineEdit()
        self.equation_edit.setPlaceholderText("e.g.  atp_c + h2o_c --> adp_c + pi_c + h_c")
        self.gpr_edit = QLineEdit()
        self.gpr_edit.setPlaceholderText("e.g.  b0001 and b0002")
        self.subsystem_edit = QLineEdit()

        self.lower = QDoubleSpinBox()
        self.lower.setRange(-1e6, 1e6)
        self.lower.setDecimals(3)
        self.lower.setValue(0.0)
        self.upper = QDoubleSpinBox()
        self.upper.setRange(-1e6, 1e6)
        self.upper.setDecimals(3)
        self.upper.setValue(1000.0)

        hint = QLabel(
            "Use '-->' for irreversible or '<=>' for reversible reactions. "
            "Metabolites that don't exist yet are created automatically."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #5f6368;")

        form = QFormLayout()
        form.addRow("Reaction ID:", self.id_edit)
        form.addRow("Name:", self.name_edit)
        form.addRow("Equation:", self.equation_edit)
        form.addRow("Lower bound:", self.lower)
        form.addRow("Upper bound:", self.upper)
        form.addRow("Gene rule (GPR):", self.gpr_edit)
        form.addRow("Subsystem:", self.subsystem_edit)

        # Live mass/charge balance indicator + one-click H+/H2O fixer.
        self.balance_label = QLabel("Balance: —")
        self.balance_label.setTextInteractionFlags(self.balance_label.textInteractionFlags())
        self.balance_btn = QPushButton("Balance with H⁺/H₂O")
        self.balance_btn.setToolTip("Close a hydrogen/oxygen/charge imbalance by adding "
                                    "water and protons (aqueous, ~pH 7). You confirm the "
                                    "exact change before it is applied.")
        self.balance_btn.clicked.connect(self._balance_clicked)
        self.balance_btn.setEnabled(False)
        # A balance that cannot be checked is a database gap, not a verdict — offer to
        # close it here rather than leaving the user with an unexplained "unknown".
        self.fetch_formula_btn = QPushButton("Fetch missing formula")
        self.fetch_formula_btn.setToolTip("Look up the missing chemical formula from "
                                          "BiGG, then from the metabolite's KEGG/ChEBI/"
                                          "InChIKey cross-references.")
        self.fetch_formula_btn.clicked.connect(self._fetch_formula_clicked)
        self.fetch_formula_btn.setVisible(False)
        balance_row = QHBoxLayout()
        balance_row.addWidget(self.balance_label, 1)
        balance_row.addWidget(self.fetch_formula_btn)
        balance_row.addWidget(self.balance_btn)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addLayout(balance_row)
        layout.addWidget(hint)
        layout.addWidget(buttons)

        if reaction is not None:
            self._populate(reaction)
        self._update_balance()

    def _populate(self, rxn: cobra.Reaction) -> None:
        self.id_edit.setText(rxn.id)
        self.id_edit.setReadOnly(True)  # don't allow renaming an existing reaction id
        self.name_edit.setText(rxn.name or "")
        self.equation_edit.setText(rxn.build_reaction_string(use_metabolite_names=False))
        self.equation_edit.setReadOnly(True)  # equation editing handled separately
        self.gpr_edit.setText(rxn.gene_reaction_rule or "")
        self.subsystem_edit.setText(rxn.subsystem or "")
        self.lower.setValue(rxn.lower_bound)
        self.upper.setValue(rxn.upper_bound)

    @property
    def is_editing(self) -> bool:
        return self._editing

    # ----- mass/charge balance -----------------------------------------
    def _update_balance(self) -> None:
        rxn = self._reaction
        self.fetch_formula_btn.setVisible(False)
        if rxn is None or getattr(rxn, "model", None) is None:
            self.balance_label.setText("Balance: — (open an existing reaction to check)")
            self.balance_btn.setEnabled(False)
            return
        verdict = balancing.assess_reaction(rxn)
        if not verdict.checkable:
            # Not "unbalanced": the missing participant's atoms were never counted, so
            # there is no verdict to give. Name the culprit — it is the fixable part.
            self.balance_label.setText(
                "Balance: cannot be checked — " + verdict.missing_formula_sentence())
            self.balance_label.setStyleSheet("color:#B06000;")
            self.balance_btn.setEnabled(False)
            self.fetch_formula_btn.setVisible(True)
            return
        if verdict.balanced:
            self.balance_label.setText("Balance: ✓ mass and charge balanced")
            self.balance_label.setStyleSheet("color:#1E8E3E; font-weight:600;")
            self.balance_btn.setEnabled(False)
            return
        self.balance_label.setText("Balance: ✗ unbalanced — "
                                   + balancing.format_residual(verdict.residual))
        self.balance_label.setStyleSheet("color:#D93025; font-weight:600;")
        self.balance_btn.setEnabled(balancing.plan_water_proton_changes(rxn) is not None)

    def _fetch_formula_clicked(self) -> None:
        """Fetch the formulas this reaction's balance check is missing, then re-check."""
        from ...core import chemistry
        from ..widgets.busy import run_busy, was_cancelled

        rxn = self._reaction
        model = getattr(rxn, "model", None) if rxn is not None else None
        if model is None:
            return
        missing = [m.id for m in balancing.missing_formula_metabolites(rxn)]
        if not missing:
            self._update_balance()
            return
        # Off the GUI thread: the lookup reaches BiGG and PubChem, and a frozen dialog
        # during a network call reads as a crash.
        ok, report = run_busy(
            self, f"Looking up {len(missing)} formula(s)…",
            lambda: chemistry.fetch_missing_formulas(model, missing, online=True),
            title="Fetch missing formula", cancelable=True)
        if not ok:
            if not was_cancelled(report):
                QMessageBox.warning(self, "Could not fetch formulas", str(report))
            return
        self._update_balance()
        QMessageBox.information(
            self, "Formula lookup",
            report.sentence() + "\n\n"
            + f"{rxn.id}: {balancing.assess_reaction(rxn).sentence()}.")

    def _balance_clicked(self) -> None:
        rxn = self._reaction
        changes = balancing.plan_water_proton_changes(rxn) if rxn is not None else None
        if not changes:
            QMessageBox.information(
                self, "Cannot auto-balance",
                "This residual cannot be closed with water and protons alone "
                "(it involves other elements or is inconsistent). It likely needs a "
                "missing substrate/cofactor or a corrected formula — fix it manually.")
            return
        summary = "\n".join(
            f"  • {'add' if c['delta'] > 0 else 'add'} {abs(c['delta']):g} "
            f"{c['metabolite_id']} to the {c['side']} side" for c in changes)
        if QMessageBox.question(
                self, "Confirm balancing changes",
                "The following will be added to balance the reaction (aqueous, ~pH 7):\n\n"
                + summary + "\n\nApply these changes?",
                QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return
        # Stage the changes and reflect them in the shown equation (applied on OK).
        self._balance_changes = changes
        self._apply_changes_to_display(changes)
        QMessageBox.information(self, "Balancing staged",
                                "The reaction will be balanced when you click OK.")

    def _apply_changes_to_display(self, changes) -> None:
        """Update the read-only equation preview to reflect staged additions."""
        # Recompute a display string by applying deltas to a copy of the stoichiometry.
        rxn = self._reaction
        coeffs = {m.id: c for m, c in rxn.metabolites.items()}
        for ch in changes:
            coeffs[ch["metabolite_id"]] = coeffs.get(ch["metabolite_id"], 0.0) + ch["delta"]
        left = " + ".join(f"{-c:g} {mid}" for mid, c in coeffs.items() if c < 0)
        right = " + ".join(f"{c:g} {mid}" for mid, c in coeffs.items() if c > 0)
        arrow = "<=>" if rxn.reversibility else "-->"
        self.equation_edit.setText(f"{left} {arrow} {right}")
        # Re-check against the staged result for user feedback.
        self.balance_label.setText("Balance: ✓ will be balanced on OK")
        self.balance_label.setStyleSheet("color:#1E8E3E; font-weight:600;")
        self.balance_btn.setEnabled(False)

    def values(self) -> dict:
        return {
            "reaction_id": self.id_edit.text().strip(),
            "name": self.name_edit.text().strip(),
            "reaction_string": self.equation_edit.text().strip(),
            "gene_reaction_rule": self.gpr_edit.text().strip(),
            "subsystem": self.subsystem_edit.text().strip(),
            "lower_bound": self.lower.value(),
            "upper_bound": self.upper.value(),
            "balance_changes": list(self._balance_changes),
        }
