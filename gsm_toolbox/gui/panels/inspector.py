"""Right-dock Inspector: details for the currently selected reaction/metabolite/gene."""

from __future__ import annotations

import cobra
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFormLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ...core import databases


class Inspector(QScrollArea):
    # Emitted when the user asks to look up enzymes for a reaction's EC numbers.
    enzymes_requested = Signal(list, str)  # [ec, …], reaction id

    def __init__(self):
        super().__init__()
        self.setWidgetResizable(True)
        self._content = QWidget()
        self._layout = QVBoxLayout(self._content)
        self._layout.setAlignment(Qt.AlignTop)
        self.setWidget(self._content)
        self.show_placeholder()

    def _clear(self) -> None:
        while self._layout.count():
            item = self._layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def show_placeholder(self) -> None:
        self._clear()
        label = QLabel("Select a reaction, metabolite or gene to see details.")
        label.setWordWrap(True)
        label.setStyleSheet("color: #5f6368; padding: 8px;")
        self._layout.addWidget(label)

    def _header(self, text: str) -> None:
        label = QLabel(text)
        label.setStyleSheet("font-weight: 700; font-size: 13px; padding: 4px 0;")
        self._layout.addWidget(label)

    def _form(self, rows: list[tuple[str, str]]) -> None:
        container = QWidget()
        form = QFormLayout(container)
        for key, value in rows:
            v = QLabel(str(value))
            v.setWordWrap(True)
            v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            form.addRow(f"{key}:", v)
        self._layout.addWidget(container)

    def show_reaction(self, rxn: cobra.Reaction, flux=None) -> None:
        self._clear()
        self._header(f"Reaction · {rxn.id}")
        ec = databases.reaction_ec_numbers(rxn)
        from ...core.network_graph import clean_label, display_reaction_name
        rows = [
            ("Name", display_reaction_name(rxn) or "—"),
            ("Equation", clean_label(rxn.build_reaction_string(use_metabolite_names=True))),
            ("Bounds", f"[{rxn.lower_bound:g}, {rxn.upper_bound:g}]"),
            ("Reversible", "yes" if rxn.reversibility else "no"),
            ("Subsystem", rxn.subsystem or "—"),
            ("EC number", ", ".join(ec) if ec else "—"),
            ("Gene rule", rxn.gene_reaction_rule or "—"),
        ]
        if flux is not None:
            rows.append(("Current flux", f"{flux:.6g}"))
        self._form(rows)
        if ec:
            btn = QPushButton(f"Find enzymes in UniProt (EC {ec[0]})")
            btn.setToolTip("List candidate enzymes and source organisms for this "
                           "reaction's EC number(s) from UniProt.")
            btn.clicked.connect(lambda _=False, e=ec, rid=rxn.id: self.enzymes_requested.emit(e, rid))
            self._layout.addWidget(btn)

    def show_metabolite(self, met: cobra.Metabolite) -> None:
        self._clear()
        self._header(f"Metabolite · {met.id}")
        self._form([
            ("Name", met.name or "—"),
            ("Formula", met.formula or "—"),
            ("Charge", str(met.charge) if met.charge is not None else "—"),
            ("Compartment", met.compartment or "—"),
            ("# Reactions", str(len(met.reactions))),
        ])

    def show_gene(self, gene: cobra.Gene) -> None:
        self._clear()
        self._header(f"Gene · {gene.id}")
        self._form([
            ("Name", gene.name or "—"),
            ("# Reactions", str(len(gene.reactions))),
            ("Reactions", ", ".join(r.id for r in gene.reactions) or "—"),
        ])
