"""The central Model Explorer: searchable Reactions / Metabolites / Genes tables.

This is the main center panel. The reactions table supports multi-selection, a
right-click context menu (show in map, edit, add to category), and programmatic
selection so the network map can sync the highlighted reaction back here.
"""

from __future__ import annotations

from typing import List, Optional

import cobra
from PySide6.QtCore import QItemSelectionModel, QSortFilterProxyModel, Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QLineEdit,
    QMenu,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..models.table_models import (
    GeneTableModel,
    MetaboliteTableModel,
    ReactionTableModel,
    SearchRole,
)


class _SearchProxy(QSortFilterProxyModel):
    """Filters rows against each item's rich ``SearchRole`` blob.

    The blob spans more than the visible columns (e.g. a reaction's gene rule,
    the genes/metabolites it touches, subsystem and annotations), so a search
    matches related terms, not just what is shown. Whitespace-separated words
    are combined with AND, so "glycolysis pyruvate" narrows progressively.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._terms: List[str] = []

    def set_query(self, text: str) -> None:
        self._terms = (text or "").lower().split()
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row: int, source_parent) -> bool:  # noqa: N802
        if not self._terms:
            return True
        index = self.sourceModel().index(source_row, 0, source_parent)
        blob = self.sourceModel().data(index, SearchRole) or ""
        return all(term in blob for term in self._terms)


class _FilterableTable(QWidget):
    """A search box above a table view that matches ids, names and related
    terms (subsystem, gene rule, neighbours, annotations)."""

    selection_changed = Signal(int)          # source-model row of current item
    context_requested = Signal(int, object)  # source row, global QPoint

    def __init__(self, source_model, placeholder: str, multi: bool = False):
        super().__init__()
        self.source_model = source_model
        self.proxy = _SearchProxy(self)
        self.proxy.setSourceModel(source_model)

        self.search = QLineEdit()
        self.search.setPlaceholderText(placeholder)
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.proxy.set_query)

        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableView.SelectRows)
        self.table.setSelectionMode(
            QAbstractItemView.ExtendedSelection if multi else QAbstractItemView.SingleSelection
        )
        self.table.setEditTriggers(QTableView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.selectionModel().selectionChanged.connect(self._emit_selection)

        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._emit_context)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.search)
        layout.addWidget(self.table)

    def _emit_selection(self, *_):
        row = self.current_source_row()
        if row is not None:
            self.selection_changed.emit(row)

    def _emit_context(self, pos):
        index = self.table.indexAt(pos)
        if not index.isValid():
            return
        src_row = self.proxy.mapToSource(index).row()
        self.context_requested.emit(src_row, self.table.viewport().mapToGlobal(pos))

    def current_source_row(self) -> Optional[int]:
        index = self.table.currentIndex()
        if not index.isValid():
            return None
        return self.proxy.mapToSource(index).row()

    def selected_source_rows(self) -> List[int]:
        return [self.proxy.mapToSource(i).row() for i in self.table.selectionModel().selectedRows()]

    def select_source_row(self, src_row: int) -> None:
        proxy_index = self.proxy.mapFromSource(self.source_model.index(src_row, 0))
        if proxy_index.isValid():
            self.table.selectionModel().setCurrentIndex(
                proxy_index,
                QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows,
            )
            self.table.scrollTo(proxy_index, QAbstractItemView.PositionAtCenter)

    def refresh_layout(self) -> None:
        self.table.resizeColumnsToContents()


class ModelExplorer(QWidget):
    """Tabbed explorer with Reactions / Metabolites / Genes (the center panel)."""

    reaction_selected = Signal(object)
    metabolite_selected = Signal(object)
    metabolite_details_requested = Signal(object)    # cobra.Metabolite -> popup + panel
    edit_metabolite_requested = Signal(object)        # cobra.Metabolite -> rename dialog
    gene_selected = Signal(object)
    reaction_activated = Signal(object)              # double-click -> edit
    show_in_map_requested = Signal(str, int)         # reaction_id, steps
    edit_reaction_requested = Signal(object)         # cobra.Reaction
    reaction_details_requested = Signal(object)      # cobra.Reaction -> info popup (incl. structures)
    add_to_category_requested = Signal(list)         # list of reaction ids

    def __init__(self):
        super().__init__()
        self.reaction_model = ReactionTableModel()
        self.metabolite_model = MetaboliteTableModel()
        self.gene_model = GeneTableModel()

        self.rxn_tab = _FilterableTable(
            self.reaction_model, "Search reactions — id, name, subsystem, gene rule…", multi=True)
        self.met_tab = _FilterableTable(
            self.metabolite_model, "Search metabolites — id, name, formula, reaction…")
        self.gene_tab = _FilterableTable(
            self.gene_model, "Search genes — id, name, reaction, subsystem…")

        self.tabs = QTabWidget()
        self.tabs.addTab(self.rxn_tab, "Reactions")
        self.tabs.addTab(self.met_tab, "Metabolites")
        self.tabs.addTab(self.gene_tab, "Genes")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.tabs)

        self.rxn_tab.selection_changed.connect(self._on_rxn_selected)
        self.met_tab.selection_changed.connect(self._on_met_selected)
        self.gene_tab.selection_changed.connect(self._on_gene_selected)
        self.rxn_tab.table.doubleClicked.connect(self._on_rxn_activated)
        self.rxn_tab.context_requested.connect(self._on_rxn_context)
        self.met_tab.context_requested.connect(self._on_met_context)
        self.gene_tab.context_requested.connect(self._on_gene_context)

    def set_model(self, model: cobra.Model) -> None:
        self.reaction_model.set_model(model)
        self.metabolite_model.set_model(model)
        self.gene_model.set_model(model)
        for tab in (self.rxn_tab, self.met_tab, self.gene_tab):
            tab.refresh_layout()
        self._update_tab_titles()

    def _update_tab_titles(self) -> None:
        self.tabs.setTabText(0, f"Reactions ({self.reaction_model.rowCount()})")
        self.tabs.setTabText(1, f"Metabolites ({self.metabolite_model.rowCount()})")
        self.tabs.setTabText(2, f"Genes ({self.gene_model.rowCount()})")

    def set_fluxes(self, fluxes: dict) -> None:
        self.reaction_model.set_fluxes(fluxes)

    def selected_reaction(self) -> Optional[cobra.Reaction]:
        row = self.rxn_tab.current_source_row()
        return None if row is None else self.reaction_model.reaction_at(row)

    def selected_reaction_ids(self) -> List[str]:
        rxns = [self.reaction_model.reaction_at(r) for r in self.rxn_tab.selected_source_rows()]
        return [r.id for r in rxns if r is not None]

    def select_reaction_by_id(self, reaction_id: str) -> None:
        """Select & scroll to a reaction (used when a node is clicked in the map)."""
        for row in range(self.reaction_model.rowCount()):
            rxn = self.reaction_model.reaction_at(row)
            if rxn is not None and rxn.id == reaction_id:
                self.tabs.setCurrentWidget(self.rxn_tab)
                self.rxn_tab.select_source_row(row)
                return

    # signal relays --------------------------------------------------
    def _on_rxn_selected(self, row: int):
        rxn = self.reaction_model.reaction_at(row)
        if rxn is not None:
            self.reaction_selected.emit(rxn)

    def _on_met_selected(self, row: int):
        met = self.metabolite_model.metabolite_at(row)
        if met is not None:
            self.metabolite_selected.emit(met)

    def _on_gene_selected(self, row: int):
        gene = self.gene_model.gene_at(row)
        if gene is not None:
            self.gene_selected.emit(gene)

    def _on_rxn_activated(self, proxy_index):
        row = self.rxn_tab.proxy.mapToSource(proxy_index).row()
        rxn = self.reaction_model.reaction_at(row)
        if rxn is not None:
            self.edit_reaction_requested.emit(rxn)

    def _on_rxn_context(self, src_row: int, global_pos):
        rxn = self.reaction_model.reaction_at(src_row)
        if rxn is None:
            return
        selected_ids = self.selected_reaction_ids() or [rxn.id]

        menu = QMenu(self)
        show_menu = menu.addMenu("Show in map")
        for steps in (1, 2, 3):
            act = show_menu.addAction(f"{steps} step{'s' if steps > 1 else ''} around")
            act.triggered.connect(lambda _=False, s=steps, rid=rxn.id: self.show_in_map_requested.emit(rid, s))
        menu.addAction("Reaction information…",
                       lambda: self.reaction_details_requested.emit(rxn))
        menu.addAction("Edit reaction…", lambda: self.edit_reaction_requested.emit(rxn))
        menu.addSeparator()
        label = f"Add {len(selected_ids)} reaction(s) to category…"
        menu.addAction(label, lambda: self.add_to_category_requested.emit(selected_ids))
        menu.exec(global_pos)

    def _on_met_context(self, src_row: int, global_pos):
        met = self.metabolite_model.metabolite_at(src_row)
        if met is None:
            return
        menu = QMenu(self)
        menu.addAction("Show details", lambda: self.metabolite_details_requested.emit(met))
        menu.addAction("Edit metabolite…", lambda: self.edit_metabolite_requested.emit(met))
        show_menu = menu.addMenu("Show in map")
        for steps in (1, 2, 3):
            act = show_menu.addAction(f"{steps} step{'s' if steps > 1 else ''} around")
            act.triggered.connect(lambda _=False, s=steps, mid=met.id: self.show_in_map_requested.emit(mid, s))
        connected = [r.id for r in met.reactions]
        menu.addSeparator()
        menu.addAction(f"Add {len(connected)} connected reaction(s) to category…",
                       lambda: self.add_to_category_requested.emit(connected))
        menu.exec(global_pos)

    def _on_gene_context(self, src_row: int, global_pos):
        gene = self.gene_model.gene_at(src_row)
        if gene is None:
            return
        menu = QMenu(self)
        menu.addAction("Show details", lambda: self.gene_selected.emit(gene))
        rxns = [r.id for r in gene.reactions]
        if rxns:
            map_menu = menu.addMenu("Show a reaction in map")
            for rid in rxns[:25]:
                act = map_menu.addAction(rid)
                act.triggered.connect(lambda _=False, r=rid: self.show_in_map_requested.emit(r, 1))
            menu.addSeparator()
            menu.addAction(f"Add this gene's {len(rxns)} reaction(s) to category…",
                           lambda: self.add_to_category_requested.emit(rxns))
        menu.exec(global_pos)
