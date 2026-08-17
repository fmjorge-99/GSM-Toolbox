"""Read-only Qt table models for reactions, metabolites and genes.

These are thin presenters over a ``cobra.Model``. Editing is done through dialogs
that call the core ``editing`` module and then ``refresh()`` the model here, so the
table models never mutate the model themselves.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import cobra
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt

from ...core import databases

# Role used by the explorer's search proxy to fetch a per-row, lower-cased
# "search blob" that spans more than just the visible columns (e.g. a reaction's
# gene rule and the names of the genes/metabolites it touches).
SearchRole = Qt.UserRole + 1


def _annotation_terms(obj) -> List[str]:
    """Flatten a cobra object's annotation values into plain strings."""
    terms: List[str] = []
    ann = getattr(obj, "annotation", None)
    if isinstance(ann, dict):
        for key, value in ann.items():
            terms.append(str(key))
            if isinstance(value, (list, tuple, set)):
                terms.extend(str(v) for v in value)
            else:
                terms.append(str(value))
    return terms


class ReactionTableModel(QAbstractTableModel):
    COLUMNS = ["ID", "Name", "Lower", "Upper", "Subsystem", "EC", "Flux"]

    def __init__(self, model: Optional[cobra.Model] = None):
        super().__init__()
        self._rxns: List[cobra.Reaction] = []
        self._search: List[str] = []
        self._fluxes: Dict[str, float] = {}
        if model is not None:
            self.set_model(model)

    def set_model(self, model: cobra.Model) -> None:
        self.beginResetModel()
        self._rxns = list(model.reactions)
        self._search = [self._build_search(r) for r in self._rxns]
        self.endResetModel()

    @staticmethod
    def _build_search(rxn: cobra.Reaction) -> str:
        parts = [rxn.id, rxn.name or "", rxn.subsystem or "",
                 rxn.gene_reaction_rule or ""]
        for g in rxn.genes:
            parts.append(g.id)
            if g.name:
                parts.append(g.name)
        for met in rxn.metabolites:
            parts.append(met.id)
            if met.name:
                parts.append(met.name)
        parts.extend(_annotation_terms(rxn))
        return " ".join(parts).lower()

    def set_fluxes(self, fluxes: Dict[str, float]) -> None:
        self._fluxes = dict(fluxes)
        if self._rxns:
            top_left = self.index(0, self.COLUMNS.index("Flux"))
            bottom_right = self.index(len(self._rxns) - 1, self.COLUMNS.index("Flux"))
            self.dataChanged.emit(top_left, bottom_right, [Qt.DisplayRole])

    def clear_fluxes(self) -> None:
        self.set_fluxes({})

    def reaction_at(self, row: int) -> Optional[cobra.Reaction]:
        if 0 <= row < len(self._rxns):
            return self._rxns[row]
        return None

    # Qt API -----------------------------------------------------------
    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._rxns)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return len(self.COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.COLUMNS[section]
        return None

    def data(self, index: QModelIndex, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        rxn = self._rxns[index.row()]
        col = index.column()
        if role == SearchRole:
            return self._search[index.row()]
        if role in (Qt.DisplayRole, Qt.ToolTipRole):
            if col == 0:
                return rxn.id
            if col == 1:
                from ...core.network_graph import clean_label
                return clean_label(rxn.name or "")
            if col == 2:
                return f"{rxn.lower_bound:g}"
            if col == 3:
                return f"{rxn.upper_bound:g}"
            if col == 4:
                return rxn.subsystem or ""
            if col == 5:
                return ", ".join(databases.reaction_ec_numbers(rxn))
            if col == 6:
                flux = self._fluxes.get(rxn.id)
                return "" if flux is None else f"{flux:.4g}"
        return None


class MetaboliteTableModel(QAbstractTableModel):
    COLUMNS = ["ID", "Name", "Formula", "Compartment"]

    def __init__(self, model: Optional[cobra.Model] = None):
        super().__init__()
        self._mets: List[cobra.Metabolite] = []
        self._search: List[str] = []
        if model is not None:
            self.set_model(model)

    def set_model(self, model: cobra.Model) -> None:
        self.beginResetModel()
        self._mets = list(model.metabolites)
        self._search = [self._build_search(m) for m in self._mets]
        self.endResetModel()

    @staticmethod
    def _build_search(met: cobra.Metabolite) -> str:
        parts = [met.id, met.name or "", met.formula or "", met.compartment or ""]
        for rxn in met.reactions:
            parts.append(rxn.id)
            if rxn.subsystem:
                parts.append(rxn.subsystem)
        parts.extend(_annotation_terms(met))
        return " ".join(parts).lower()

    def metabolite_at(self, row: int) -> Optional[cobra.Metabolite]:
        if 0 <= row < len(self._mets):
            return self._mets[row]
        return None

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._mets)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return len(self.COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.COLUMNS[section]
        return None

    def data(self, index: QModelIndex, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        if role == SearchRole:
            return self._search[index.row()]
        if role not in (Qt.DisplayRole, Qt.ToolTipRole):
            return None
        met = self._mets[index.row()]
        from ...core.network_graph import short_metabolite_name
        display_name = short_metabolite_name(met.id, met.name or "") if met.name else ""
        return [met.id, display_name, met.formula or "", met.compartment or ""][index.column()]


class GeneTableModel(QAbstractTableModel):
    COLUMNS = ["ID", "Name", "# Reactions"]

    def __init__(self, model: Optional[cobra.Model] = None):
        super().__init__()
        self._genes: List[cobra.Gene] = []
        self._search: List[str] = []
        if model is not None:
            self.set_model(model)

    def set_model(self, model: cobra.Model) -> None:
        self.beginResetModel()
        self._genes = list(model.genes)
        self._search = [self._build_search(g) for g in self._genes]
        self.endResetModel()

    @staticmethod
    def _build_search(gene: cobra.Gene) -> str:
        parts = [gene.id, gene.name or ""]
        for rxn in gene.reactions:
            parts.append(rxn.id)
            if rxn.subsystem:
                parts.append(rxn.subsystem)
        parts.extend(_annotation_terms(gene))
        return " ".join(parts).lower()

    def gene_at(self, row: int) -> Optional[cobra.Gene]:
        if 0 <= row < len(self._genes):
            return self._genes[row]
        return None

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._genes)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return len(self.COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.COLUMNS[section]
        return None

    def data(self, index: QModelIndex, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        if role == SearchRole:
            return self._search[index.row()]
        if role not in (Qt.DisplayRole, Qt.ToolTipRole):
            return None
        gene = self._genes[index.row()]
        return [gene.id, gene.name or "", str(len(gene.reactions))][index.column()]
