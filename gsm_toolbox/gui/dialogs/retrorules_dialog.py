"""Show rule-based retrosynthesis (RetroRules) suggestions for a target.

The steps are PREDICTIONS from generalised reaction rules, not database-backed
reactions, so the dialog says so plainly and renders each step as
``product ← precursors`` — with the 2-D structure of every compound drawn on the right,
and a human-readable name where one can be resolved.

Usability (refinement 3.2):
* the target's own fetched structure + name are shown at the top;
* several ranked alternative routes can be offered — a selector switches between them;
* each step has a checkbox, so the user adds only the steps they want (not all-or-nothing).
"Add selected reactions" hands the chosen steps to the main Pathway Design panel.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

_IMG = 116          # structure thumbnail size (px)


def _ec_links(ec_field: str) -> str:
    """Render a rule's EC field as clickable ExPASy ENZYME links.

    RetroRules stores one or several EC numbers per rule; each becomes a link so the user
    can jump from "this step is proposed" to "this is the enzyme class that does it" —
    the first concrete step toward finding a gene to clone.
    """
    import re
    ecs = re.findall(r"\d+\.\d+\.\d+\.(?:\d+|-)", str(ec_field or ""))
    if not ecs:
        return f"EC {ec_field}"
    seen, links = set(), []
    for ec in ecs:
        if ec in seen:
            continue
        seen.add(ec)
        links.append(f"<a href='https://enzyme.expasy.org/EC/{ec}'>{ec}</a>")
    return "EC " + ", ".join(links)


class RetroRulesOptionsDialog(QDialog):
    """Ask how to run RetroRules before searching: how many alternative routes, how deep,
    and how to rank them (refinement 3.2 — expanded input options)."""

    def __init__(self, parent, target_name: str):
        super().__init__(parent)
        from PySide6.QtWidgets import QDialogButtonBox, QFormLayout, QSpinBox
        self.setWindowTitle(f"RetroRules options — {target_name}")
        form = QFormLayout(self)

        self.n_alt = QSpinBox()
        self.n_alt.setRange(1, 8)
        self.n_alt.setValue(3)
        self.n_alt.setToolTip("How many genuinely different routes to look for. Each uses "
                              "a different first disconnection of the target.")
        form.addRow("Alternative routes:", self.n_alt)

        self.max_steps = QSpinBox()
        self.max_steps.setRange(1, 8)
        self.max_steps.setValue(4)
        self.max_steps.setToolTip("Maximum heterologous steps back from the target.")
        form.addRow("Max steps:", self.max_steps)

        self.rank = QComboBox()
        # (key, human label) — keys match retrorules.RANK_KEYS.
        self._rank_keys = ["score", "steps", "precursors"]
        self.rank.addItems([
            "RetroRules reliability score (recommended)",
            "Fewest steps",
            "Fewest native precursors",
        ])
        self.rank.setToolTip("How to order the alternative routes that are found.")
        form.addRow("Rank routes by:", self.rank)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        form.addRow(bb)

    def values(self):
        return {
            "n_alternatives": self.n_alt.value(),
            "max_steps": self.max_steps.value(),
            "rank_by": self._rank_keys[self.rank.currentIndex()],
        }


def _structure_pixmap(smiles: str, cache: Dict[str, QPixmap],
                      size: int = _IMG) -> Optional[QPixmap]:
    ck = f"{smiles}@{size}"
    if ck in cache:
        return cache[ck]
    pm: Optional[QPixmap] = None
    try:
        from ..widgets.structure_fetcher import _render_bw
        png = _render_bw(smiles=smiles, size=size)
        if png:
            p = QPixmap()
            if p.loadFromData(png):
                pm = p
    except Exception:  # noqa: BLE001 — a missing image must never break the dialog
        pm = None
    cache[ck] = pm
    return pm


class RetroRulesDialog(QDialog):
    add_reactions_requested = Signal()      # user chose "Add selected reactions"

    def __init__(self, parent, target_name: str, routes, *, target_smiles: str = "",
                 rank_label: str = "", ruleset_note: str = "",
                 name_for: Optional[Callable[[str], str]] = None,
                 seed: Optional[int] = None):
        super().__init__(parent)
        self.setWindowTitle(f"Rule-based suggestions — {target_name}")
        self.resize(940, 680)
        from ..widgets.dialog_util import clamp_to_screen
        clamp_to_screen(self)

        # Accept either a single route (back-compat) or a ranked list of routes.
        if routes is None:
            self._routes: List = []
        elif isinstance(routes, list):
            self._routes = [r for r in routes if r is not None]
        else:
            self._routes = [routes]
        self._target_name = target_name
        self._target_smiles = target_smiles
        self._name_for = name_for or (lambda s: "")
        self._img_cache: Dict[str, QPixmap] = {}
        self._step_checks: List[QCheckBox] = []
        self._current_idx = 0

        layout = QVBoxLayout(self)

        # -- target header: its own structure + name, so the user sees WHAT is being made.
        layout.addWidget(self._build_target_header())

        warn = QLabel(
            "⚠ These steps are <b>predictions</b> from generalised reaction rules, not "
            "reactions from a database. Treat them as hypotheses: check each is real "
            "chemistry and that an enzyme exists before relying on it."
            + (f"<br><span style='color:#5f6368'>{ruleset_note}</span>" if ruleset_note
               else ""))
        warn.setWordWrap(True)
        warn.setStyleSheet("color:#8a6d00; background:#fff8e1; padding:6px; "
                           "border-radius:4px;")
        layout.addWidget(warn)

        # Reproducibility / coverage caveat (L1). The search samples a very large rule
        # set, so what comes back depends on the exploration order — which is fixed by
        # the seed. Users must know these are A sample, not THE answer.
        # Report the seed this SEARCH actually used, not whatever is currently saved in
        # preferences — otherwise changing the preference (or running a search with an
        # explicit seed) makes the dialog quote a seed that would not reproduce it, which
        # defeats the purpose of showing it at all.
        from ...core import preferences
        shown_seed = preferences.retrorules_seed() if seed is None else seed
        repro = QLabel(
            f"<b>Seed {shown_seed}</b> — re-running returns these same routes. This is "
            "a sample of possible chemistry, not an exhaustive enumeration; another "
            "seed may surface different, equally valid routes.")
        repro.setWordWrap(True)
        repro.setTextFormat(Qt.RichText)
        repro.setStyleSheet("color:#1967D2; background:#E8F0FE; padding:6px; "
                            "border-radius:4px;")
        layout.addWidget(repro)

        # -- route selector (only when there is more than one alternative).
        self._head = QLabel("")
        self._head.setWordWrap(True)
        layout.addWidget(self._head)
        if len(self._routes) > 1:
            sel_row = QHBoxLayout()
            sel_row.addWidget(QLabel("Alternative route:"))
            self.route_combo = QComboBox()
            for i, r in enumerate(self._routes, 1):
                tag = "complete" if getattr(r, "complete", False) else "partial"
                self.route_combo.addItem(
                    f"Route {i} — {r.n_steps} step(s) · score {r.score:.2f} · {tag}")
            self.route_combo.currentIndexChanged.connect(self._on_route_changed)
            sel_row.addWidget(self.route_combo, 1)
            if rank_label:
                lbl = QLabel(f"(ranked by {rank_label})")
                lbl.setStyleSheet("color:#5f6368;")
                sel_row.addWidget(lbl)
            layout.addLayout(sel_row)

        # -- scrollable body of steps (repopulated when the route changes).
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        self._body = QVBoxLayout(inner)
        self._body.setAlignment(Qt.AlignTop)
        scroll.setWidget(inner)
        layout.addWidget(scroll, 1)

        btn_row = QHBoxLayout()
        self._check_all = QPushButton("Select all steps")
        self._check_all.clicked.connect(lambda: self._set_all_checked(True))
        self._check_none = QPushButton("Clear")
        self._check_none.clicked.connect(lambda: self._set_all_checked(False))
        btn_row.addWidget(self._check_all)
        btn_row.addWidget(self._check_none)
        btn_row.addStretch(1)
        self.add_btn = QPushButton("Add selected reactions")
        self.add_btn.setObjectName("primary")
        self.add_btn.setToolTip("Load the checked steps into Pathway Design as a "
                                "suggested pathway, then use “Apply pathway” to add them "
                                "to your model (naming the reactions and metabolites).")
        self.add_btn.clicked.connect(self._on_add)
        btn_row.addWidget(self.add_btn)
        close = QPushButton("Close")
        close.clicked.connect(self.reject)
        btn_row.addWidget(close)
        layout.addLayout(btn_row)

        self._populate()

    # -- public accessors used by the main window ------------------------------
    def selected_route(self):
        return self._routes[self._current_idx] if self._routes else None

    def selected_step_indices(self) -> List[int]:
        return [i for i, cb in enumerate(self._step_checks) if cb.isChecked()]

    # -- header ----------------------------------------------------------------
    def _build_target_header(self) -> QWidget:
        w = QFrame()
        w.setFrameShape(QFrame.StyledPanel)
        h = QHBoxLayout(w)
        img = QLabel()
        img.setAlignment(Qt.AlignCenter)
        img.setFixedSize(_IMG, _IMG)
        pm = _structure_pixmap(self._target_smiles, self._img_cache) \
            if self._target_smiles else None
        if pm is not None:
            img.setPixmap(pm)
        else:
            img.setText("(no structure)")
            img.setStyleSheet("color:#9aa0a6; font-size:10px;")
        h.addWidget(img, 0)
        txt = QLabel(f"<b>Target: {self._target_name}</b>"
                     + (f"<br><span style='color:#5f6368'>{self._target_smiles}</span>"
                        if self._target_smiles else ""))
        txt.setWordWrap(True)
        txt.setTextFormat(Qt.RichText)
        h.addWidget(txt, 1)
        return w

    # -- body ------------------------------------------------------------------
    def _label_for(self, smiles: str) -> str:
        name = self._name_for(smiles)
        return f"{name} ({smiles})" if name else smiles

    def _compound_widget(self, smiles: str) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(2, 2, 2, 2)
        v.setSpacing(1)
        pm = _structure_pixmap(smiles, self._img_cache)
        img = QLabel()
        img.setAlignment(Qt.AlignCenter)
        if pm is not None:
            img.setPixmap(pm)
        else:
            img.setText("(no structure)")
            img.setStyleSheet("color:#9aa0a6; font-size:10px;")
        img.setFixedSize(_IMG, _IMG)
        v.addWidget(img, 0, Qt.AlignCenter)
        cap = QLabel(self._name_for(smiles) or smiles)
        cap.setWordWrap(True)
        cap.setAlignment(Qt.AlignCenter)
        cap.setStyleSheet("font-size:10px; color:#3c4043;")
        cap.setFixedWidth(_IMG + 20)
        v.addWidget(cap, 0, Qt.AlignCenter)
        return w

    def _clear_body(self) -> None:
        while self._body.count():
            item = self._body.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
        self._step_checks = []

    def _populate(self) -> None:
        self._clear_body()
        route = self.selected_route()
        if route is None or not route.steps:
            self._head.setText("")
            self._body.addWidget(QLabel("(no steps to show)"))
            self.add_btn.setEnabled(False)
            self._check_all.setEnabled(False)
            self._check_none.setEnabled(False)
            return
        complete = bool(getattr(route, "complete", False))
        n = route.n_steps
        if complete:
            head = (f"<b>Rule-based route in {n} step(s)</b> (score {route.score:.2f}), "
                    "grounded in compounds your host already makes. Tick the steps to add.")
        else:
            head = (f"<b>Partial route ({n} step(s), score {route.score:.2f}).</b> Some "
                    "branches did not reach a native precursor within the budget. Tick "
                    "the steps to add.")
        self._head.setText(head)

        for i, s in enumerate(route.steps):
            row = QFrame()
            row.setFrameShape(QFrame.StyledPanel)
            h = QHBoxLayout(row)
            cb = QCheckBox()
            cb.setChecked(True)
            cb.setToolTip("Include this step when adding the pathway.")
            self._step_checks.append(cb)
            h.addWidget(cb, 0, Qt.AlignTop)
            left = QLabel(
                f"<b>{i + 1}. {self._label_for(s['product'])}</b>"
                f"<br>&nbsp;&nbsp;← {' + '.join(self._label_for(p) for p in s['precursors'])}"
                f"<br><span style='color:#5f6368'>rule {s.get('rule_id','')}"
                + (f" · {_ec_links(s['ec'])}" if s.get('ec') else "")
                + (f" · score {s.get('score', 0.0):.2f}" if s.get('score') else "")
                + "</span>")
            left.setWordWrap(True)
            left.setTextFormat(Qt.RichText)
            # EC numbers become links to ExPASy ENZYME, so the user can go straight from
            # a proposed step to the enzyme class that performs it.
            left.setOpenExternalLinks(True)
            h.addWidget(left, 1)
            imgs = QHBoxLayout()
            for p in s["precursors"]:
                imgs.addWidget(self._compound_widget(p))
            arrow = QLabel("→")
            arrow.setStyleSheet("font-size:20px; color:#5f6368;")
            imgs.addWidget(arrow, 0, Qt.AlignCenter)
            imgs.addWidget(self._compound_widget(s["product"]))
            imgs_w = QWidget()
            imgs_w.setLayout(imgs)
            h.addWidget(imgs_w, 0)
            self._body.addWidget(row)

        if route.terminal_precursors:
            self._body.addWidget(QLabel("<b>Native precursors the route starts from:</b>"))
            nat = QHBoxLayout()
            for p in route.terminal_precursors:
                nat.addWidget(self._compound_widget(p))
            nat.addStretch(1)
            nat_w = QWidget()
            nat_w.setLayout(nat)
            self._body.addWidget(nat_w)

        self.add_btn.setEnabled(True)
        self._check_all.setEnabled(True)
        self._check_none.setEnabled(True)

    def _on_route_changed(self, idx: int) -> None:
        self._current_idx = idx
        self._populate()

    def _set_all_checked(self, on: bool) -> None:
        for cb in self._step_checks:
            cb.setChecked(on)

    def _on_add(self) -> None:
        if not self.selected_step_indices():
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.information(self, "No steps selected",
                                    "Tick at least one step to add.")
            return
        self.add_reactions_requested.emit()
        self.accept()
