"""Write and edit regulatory rules against the loaded model.

A rule is four decisions: *what is sensed*, *how sharply the cell responds*, *what is
affected*, and *how much*. This dialog asks exactly those four and nothing else, because
the JSON behind them is easy to get subtly wrong by hand — a misspelled reaction id
produces a rule that loads, validates, fires, and changes nothing.

Two things the editor does that a text editor cannot:

* **targets are picked from the model**, so a rule cannot name a reaction that does not
  exist, and a whole pathway can be selected by subsystem in one action;
* **the rule is described back in prose** as it is edited ("fires fully below 1 mM"),
  which is the fastest way to catch a response curve pointing the wrong way — the single
  most common mistake, and one that produces plausible-looking output.

Confidence is a required field rather than an optional one. A rule set whose thresholds
are guesses is still useful; a rule set that does not say which thresholds are guesses is
not, because every number in it reads as measured.
"""
from __future__ import annotations

import copy
from typing import List, Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QGroupBox, QHBoxLayout, QInputDialog, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMessageBox, QPushButton, QScrollArea, QSplitter, QTextEdit,
    QVBoxLayout, QWidget)

from ...core import regulation as reg
from ...core import rule_library as lib
from .. import style

#: Effects, with the plain-language description shown next to each.
EFFECTS = [
    (reg.GATE, "Switch reactions on or off"),
    (reg.SCALE, "Scale reaction capacity"),
    (reg.ENZYME_COST, "Change enzyme cost (enzyme-constrained models)"),
    (reg.BUDGET, "Change the total protein budget"),
    (reg.BIOMASS, "Change a biomass component"),
    (reg.PARAMETER, "Change a simulation parameter"),
]

RESPONSE_KINDS = [
    ("hill", "Saturating (Hill) — the usual shape for an induced or repressed system"),
    ("step", "Hard switch — only where the biology really is all-or-nothing"),
    ("ramp", "Linear ramp — honest when the shape is unknown"),
]

CONFIDENCES = [
    (reg.MEASURED, "measured — threshold and shape from primary literature"),
    (reg.INFERRED, "inferred — mechanism established, quantities from related work"),
    (reg.ASSUMED, "assumed — mechanism established, quantities chosen for plausibility"),
]


class RuleEditorDialog(QDialog):
    """Create or edit a rule set for the loaded model."""

    def __init__(self, parent, ruleset: reg.RuleSet, model=None, path: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Regulatory rule editor")
        from ..widgets.dialog_util import clamp_to_screen

        # Size against the window this was opened from, not a fixed number. The editor
        # stacks four groups of fields, and on a smaller display the fixed 1080x660 was
        # taller than the app itself — the lower fields were simply cut off, which reads
        # as "the panels do not respond" because the controls are not reachable.
        available = None
        top = parent.window() if parent is not None else None
        if top is not None:
            available = top.size()
        screen = self.screen() or QGuiApplication.primaryScreen()
        limit = screen.availableGeometry() if screen is not None else None
        width, height = 1080, 660
        if available is not None:
            width = min(width, max(720, available.width() - 80))
            height = min(height, max(480, available.height() - 80))
        if limit is not None:
            width = min(width, limit.width() - 40)
            height = min(height, limit.height() - 60)
        self.resize(width, height)

        self._ruleset = copy.deepcopy(ruleset)
        self._model = model
        self._path = path
        self._loading = False
        self.saved_path = ""

        outer = QVBoxLayout(self)
        outer.addWidget(self._header())

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._rule_list())

        # Scrollable, so every field stays reachable however short the dialog is. Without
        # this the lower groups were clipped off the bottom and could not be edited.
        editor_scroll = QScrollArea()
        editor_scroll.setWidget(self._editor())
        editor_scroll.setWidgetResizable(True)
        editor_scroll.setFrameShape(QScrollArea.NoFrame)
        editor_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        splitter.addWidget(editor_scroll)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, max(420, self.width() - 300)])
        outer.addWidget(splitter, 1)

        self.fit_note = QLabel("")
        self.fit_note.setWordWrap(True)
        self.fit_note.setTextFormat(Qt.RichText)
        outer.addWidget(self.fit_note)

        buttons = QDialogButtonBox()
        buttons.addButton("Save to library", QDialogButtonBox.AcceptRole
                          ).clicked.connect(self._save)
        buttons.addButton("Save as…", QDialogButtonBox.ActionRole
                          ).clicked.connect(self._save_as)
        buttons.addButton(QDialogButtonBox.Close).clicked.connect(self.reject)
        outer.addWidget(buttons)

        self._reload_list()
        self._check_fit()
        clamp_to_screen(self)

    # -- layout ---------------------------------------------------------------------
    def _header(self) -> QWidget:
        box = QGroupBox("Rule set")
        form = QFormLayout(box)
        self.set_name = QLineEdit(self._ruleset.name)
        form.addRow("Name:", self.set_name)
        self.set_organism = QLineEdit(getattr(self._ruleset, "organism", ""))
        self.set_organism.setPlaceholderText("Organism this set describes")
        self.set_organism.setToolTip(
            "Regulation is organism-specific. Recording the organism is what stops this "
            "set being applied to a model it does not describe.")
        form.addRow("Organism:", self.set_organism)
        self.set_description = QLineEdit(getattr(self._ruleset, "description", ""))
        form.addRow("Description:", self.set_description)
        return box

    def _rule_list(self) -> QWidget:
        panel = QWidget()
        v = QVBoxLayout(panel)
        v.setContentsMargins(0, 0, 0, 0)
        v.addWidget(QLabel("<b>Rules</b>"))
        self.rules = QListWidget()
        self.rules.currentRowChanged.connect(self._select_rule)
        self.rules.itemChanged.connect(self._toggle_enabled)
        v.addWidget(self.rules, 1)

        row = QHBoxLayout()
        for label, slot, tip in (
                ("Add", self._add_rule, "Start a new rule."),
                ("Duplicate", self._duplicate_rule,
                 "Copy the selected rule — the quick way to write a second threshold for "
                 "the same regulator."),
                ("Delete", self._delete_rule, "Remove the selected rule.")):
            button = QPushButton(label)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            row.addWidget(button)
        v.addLayout(row)
        return panel

    def _editor(self) -> QWidget:
        panel = QWidget()
        v = QVBoxLayout(panel)
        v.setContentsMargins(0, 0, 0, 0)

        what = QGroupBox("What the rule senses")
        form = QFormLayout(what)
        self.name = QLineEdit()
        self.name.editingFinished.connect(self._commit)
        form.addRow("Rule name:", self.name)
        self.regulator = QLineEdit()
        self.regulator.setPlaceholderText("Regulator name")
        self.regulator.setToolTip("The protein or system responsible, for the record.")
        self.regulator.editingFinished.connect(self._commit)
        form.addRow("Regulator:", self.regulator)

        self.sensor = QComboBox()
        self.sensor.setEditable(True)
        for key, label in reg.STANDARD_SENSORS.items():
            self.sensor.addItem(label, key)
        self.sensor.setToolTip(
            "What the rule reads. Named sensors survive a change of nutrient "
            "source, so a rule stays valid across sources of the same element.")
        self.sensor.currentTextChanged.connect(lambda _: self._commit())
        form.addRow("Sensor:", self.sensor)
        v.addWidget(what)

        how = QGroupBox("How it responds")
        hform = QFormLayout(how)
        self.kind = QComboBox()
        for key, label in RESPONSE_KINDS:
            self.kind.addItem(label, key)
        self.kind.currentIndexChanged.connect(self._kind_changed)
        hform.addRow("Shape:", self.kind)

        self.param_a = QDoubleSpinBox()
        self.param_b = QDoubleSpinBox()
        for spin in (self.param_a, self.param_b):
            spin.setDecimals(4)
            spin.setRange(-1e6, 1e6)
            spin.setMaximumWidth(140)
            spin.valueChanged.connect(lambda _: self._commit())
        self.label_a = QLabel("Threshold:")
        self.label_b = QLabel("")
        hform.addRow(self.label_a, self.param_a)
        hform.addRow(self.label_b, self.param_b)

        self.rising = QCheckBox("Fires as the sensor rises")
        self.rising.setToolTip(
            "Uncheck for de-repression — a system that switches ON as the sensor FALLS, "
            "which is how most nutrient-scarcity regulators behave.")
        self.rising.toggled.connect(lambda _: self._commit())
        hform.addRow("", self.rising)

        self.reads = QLabel("")
        self.reads.setWordWrap(True)
        self.reads.setStyleSheet(f"color:{style.TEXT_MUTED};")
        hform.addRow("In words:", self.reads)
        v.addWidget(how)

        does = QGroupBox("What it affects")
        dform = QFormLayout(does)
        self.effect = QComboBox()
        for key, label in EFFECTS:
            self.effect.addItem(label, key)
        self.effect.currentIndexChanged.connect(self._effect_changed)
        dform.addRow("Effect:", self.effect)

        target_row = QHBoxLayout()
        self.targets = QLineEdit()
        self.targets.setPlaceholderText("no targets selected")
        self.targets.editingFinished.connect(self._commit)
        target_row.addWidget(self.targets, 1)
        self.pick_btn = QPushButton("Choose…")
        self.pick_btn.setToolTip(
            "Pick from the loaded model — by name, or a whole subsystem at once.")
        self.pick_btn.clicked.connect(self._pick_targets)
        target_row.addWidget(self.pick_btn)
        dform.addRow("Targets:", target_row)

        self.magnitude = QDoubleSpinBox()
        self.magnitude.setDecimals(4)
        self.magnitude.setRange(0.0, 1000.0)
        self.magnitude.setMaximumWidth(140)
        self.magnitude.valueChanged.connect(lambda _: self._commit())
        self.magnitude_hint = QLabel("")
        self.magnitude_hint.setWordWrap(True)
        self.magnitude_hint.setStyleSheet(f"color:{style.TEXT_MUTED};")
        mrow = QHBoxLayout()
        mrow.addWidget(self.magnitude)
        mrow.addWidget(self.magnitude_hint, 1)
        dform.addRow("Magnitude:", mrow)
        v.addWidget(does)

        evidence = QGroupBox("Evidence")
        eform = QFormLayout(evidence)
        self.confidence = QComboBox()
        for key, label in CONFIDENCES:
            self.confidence.addItem(label, key)
        self.confidence.setToolTip(
            "How well founded the numbers are. Results that depend on an 'assumed' rule "
            "are flagged wherever they appear.")
        self.confidence.currentIndexChanged.connect(lambda _: self._commit())
        eform.addRow("Confidence:", self.confidence)
        self.provenance = QLineEdit()
        self.provenance.setPlaceholderText("citation, or how the threshold was chosen")
        self.provenance.editingFinished.connect(self._commit)
        eform.addRow("Provenance:", self.provenance)
        self.description = QTextEdit()
        self.description.setMaximumHeight(60)
        self.description.textChanged.connect(self._commit)
        eform.addRow("Description:", self.description)
        v.addWidget(evidence)
        v.addStretch(1)

        self._editor_widgets = panel
        panel.setEnabled(False)
        return panel

    # -- rule list ------------------------------------------------------------------
    def _reload_list(self, keep: int = -1) -> None:
        self._loading = True
        self.rules.clear()
        for rule in self._ruleset.rules:
            item = QListWidgetItem(rule.name or "(unnamed)")
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if rule.enabled else Qt.Unchecked)
            item.setToolTip(f"{rule.regulator} · {rule.effect} · {rule.confidence}")
            self.rules.addItem(item)
        self._loading = False
        if self._ruleset.rules:
            self.rules.setCurrentRow(min(max(keep, 0), len(self._ruleset.rules) - 1))
        else:
            self._editor_widgets.setEnabled(False)

    def _current(self) -> Optional[reg.Rule]:
        row = self.rules.currentRow()
        if 0 <= row < len(self._ruleset.rules):
            return self._ruleset.rules[row]
        return None

    def _toggle_enabled(self, item: QListWidgetItem) -> None:
        if self._loading:
            return
        row = self.rules.row(item)
        if 0 <= row < len(self._ruleset.rules):
            self._ruleset.rules[row].enabled = item.checkState() == Qt.Checked

    def _add_rule(self) -> None:
        name, ok = QInputDialog.getText(self, "New rule", "Rule name:")
        if not ok or not name.strip():
            return
        spec = {"kind": "hill", "half": 1.0, "coefficient": 2.0, "rising": True}
        rule = reg.Rule(
            name=name.strip(), regulator="", sensor="nitrogen_mM",
            response=reg.build_response(spec), effect=reg.GATE, targets=[],
            magnitude=0.0, confidence=reg.ASSUMED, spec=spec)
        self._ruleset.rules.append(rule)
        self._reload_list(keep=len(self._ruleset.rules) - 1)
        self._check_fit()

    def _duplicate_rule(self) -> None:
        rule = self._current()
        if rule is None:
            return
        clone = copy.deepcopy(rule)
        clone.name = f"{rule.name} (copy)"
        self._ruleset.rules.insert(self.rules.currentRow() + 1, clone)
        self._reload_list(keep=self.rules.currentRow() + 1)

    def _delete_rule(self) -> None:
        row = self.rules.currentRow()
        if not (0 <= row < len(self._ruleset.rules)):
            return
        name = self._ruleset.rules[row].name
        if QMessageBox.question(self, "Delete rule", f"Delete '{name}'?") \
                != QMessageBox.Yes:
            return
        del self._ruleset.rules[row]
        self._reload_list(keep=row - 1)
        self._check_fit()

    # -- editing --------------------------------------------------------------------
    def _select_rule(self, row: int) -> None:
        rule = self._current()
        self._editor_widgets.setEnabled(rule is not None)
        if rule is None:
            return
        self._loading = True
        self.name.setText(rule.name)
        self.regulator.setText(rule.regulator)

        index = self.sensor.findData(rule.sensor)
        if index >= 0:
            self.sensor.setCurrentIndex(index)
        else:
            self.sensor.setEditText(rule.sensor)

        spec = rule.spec or {"kind": "step", "threshold": 0.0}
        kind_index = self.kind.findData(spec.get("kind", "step"))
        self.kind.setCurrentIndex(max(0, kind_index))
        self._sync_parameter_labels(spec.get("kind", "step"))
        if spec.get("kind") == "hill":
            self.param_a.setValue(float(spec.get("half", 0.0)))
            self.param_b.setValue(float(spec.get("coefficient", 2.0)))
        elif spec.get("kind") == "ramp":
            self.param_a.setValue(float(spec.get("low", 0.0)))
            self.param_b.setValue(float(spec.get("high", 0.0)))
        else:
            self.param_a.setValue(float(spec.get("threshold", 0.0)))
        self.rising.setChecked(bool(spec.get("rising", spec.get("above", True))))

        self.effect.setCurrentIndex(max(0, self.effect.findData(rule.effect)))
        self.targets.setText(", ".join(rule.targets))
        self.magnitude.setValue(float(rule.magnitude))
        self.confidence.setCurrentIndex(
            max(0, self.confidence.findData(rule.confidence)))
        self.provenance.setText(rule.provenance)
        self.description.setPlainText(rule.description)
        self._loading = False
        self._sync_effect_hints()
        self._describe()

    def _sync_parameter_labels(self, kind: str) -> None:
        if kind == "hill":
            self.label_a.setText("Half-maximal at:")
            self.label_b.setText("Steepness (n):")
            self.param_b.setVisible(True)
            self.label_b.setVisible(True)
        elif kind == "ramp":
            self.label_a.setText("Starts at:")
            self.label_b.setText("Complete at:")
            self.param_b.setVisible(True)
            self.label_b.setVisible(True)
        else:
            self.label_a.setText("Threshold:")
            self.label_b.setText("")
            self.param_b.setVisible(False)
            self.label_b.setVisible(False)

    def _kind_changed(self) -> None:
        self._sync_parameter_labels(self.kind.currentData())
        self._commit()

    def _effect_changed(self) -> None:
        self._sync_effect_hints()
        self._commit()

    def _sync_effect_hints(self) -> None:
        effect = self.effect.currentData()
        hints = {
            reg.GATE: "Ignored — a gate is on or off.",
            reg.SCALE: "Capacity multiplier at full activation (0.5 halves it).",
            reg.ENZYME_COST: "Enzyme-cost multiplier at full activation.",
            reg.BUDGET: "Protein-budget multiplier at full activation.",
            reg.BIOMASS: "Biomass-coefficient multiplier at full activation.",
            reg.PARAMETER: "Parameter multiplier at full activation.",
        }
        self.magnitude_hint.setText(hints.get(effect, ""))
        self.magnitude.setEnabled(effect != reg.GATE)
        needs_model_targets = effect not in (reg.BUDGET,)
        self.pick_btn.setEnabled(needs_model_targets and self._model is not None)
        self.targets.setEnabled(needs_model_targets)
        if effect == reg.BUDGET:
            self.targets.setPlaceholderText("not needed — this acts on the whole budget")
        elif effect == reg.PARAMETER:
            self.targets.setPlaceholderText("parameter name")
        elif effect == reg.BIOMASS:
            self.targets.setPlaceholderText("metabolite ids in the biomass reaction")
        else:
            self.targets.setPlaceholderText("reaction ids — use Choose…")

    def _spec_from_widgets(self) -> dict:
        kind = self.kind.currentData()
        if kind == "hill":
            return {"kind": "hill", "half": self.param_a.value(),
                    "coefficient": self.param_b.value() or 2.0,
                    "rising": self.rising.isChecked()}
        if kind == "ramp":
            return {"kind": "ramp", "low": self.param_a.value(),
                    "high": self.param_b.value(), "rising": self.rising.isChecked()}
        return {"kind": "step", "threshold": self.param_a.value(),
                "above": self.rising.isChecked()}

    def _commit(self) -> None:
        """Write the widgets back into the selected rule."""
        if self._loading:
            return
        rule = self._current()
        if rule is None:
            return
        rule.name = self.name.text().strip() or rule.name
        rule.regulator = self.regulator.text().strip()
        rule.sensor = (self.sensor.currentData()
                       if self.sensor.findText(self.sensor.currentText()) >= 0
                       else self.sensor.currentText().strip())
        rule.spec = self._spec_from_widgets()
        try:
            rule.response = reg.build_response(rule.spec)
        except (TypeError, ValueError):
            pass                      # keep the previous callable until the spec is sane
        rule.effect = self.effect.currentData()
        rule.targets = [t.strip() for t in self.targets.text().split(",") if t.strip()]
        rule.magnitude = self.magnitude.value()
        rule.confidence = self.confidence.currentData()
        rule.provenance = self.provenance.text().strip()
        rule.description = self.description.toPlainText().strip()

        row = self.rules.currentRow()
        if 0 <= row < self.rules.count():
            self.rules.item(row).setText(rule.name or "(unnamed)")
        self._describe()

    def _describe(self) -> None:
        rule = self._current()
        if rule is None:
            self.reads.setText("")
            return
        sensor = reg.sensor_label(rule.sensor)
        self.reads.setText(f"Reading <b>{sensor}</b>, this rule "
                           f"{rule.describe_response()}.")

    # -- targets --------------------------------------------------------------------
    def _pick_targets(self) -> None:
        rule = self._current()
        if rule is None or self._model is None:
            return
        from .reaction_browser_dialog import ReactionBrowserDialog

        metabolites = rule.effect == reg.BIOMASS
        prompt = ("Choose the biomass components this rule changes."
                  if metabolites else
                  "Choose the reactions this rule acts on. Filter by subsystem and use "
                  "<b>Add all shown</b> to target a whole pathway at once.")
        picked = ReactionBrowserDialog.pick(
            self, self._model, rule.targets,
            title=f"Targets for '{rule.name}'", prompt=prompt,
            allow_metabolites=metabolites)
        if picked is None:
            return
        self.targets.setText(", ".join(picked))
        self._commit()
        self._check_fit()

    def _check_fit(self) -> None:
        """Say whether these rules can actually affect the loaded model."""
        if self._model is None:
            self.fit_note.setText(
                f"<span style='color:{style.TEXT_MUTED}'>No model loaded — targets "
                "cannot be checked against anything.</span>")
            return
        report = lib.fit(self._ruleset, self._model)
        text = report.summary()
        colour = "#c5221f" if report.inapplicable else style.TEXT_MUTED
        self.fit_note.setText(f"<span style='color:{colour}'>{text}</span>")

    # -- saving ---------------------------------------------------------------------
    def _harvest(self) -> reg.RuleSet:
        self._commit()
        self._ruleset.name = self.set_name.text().strip() or "Rule set"
        self._ruleset.organism = self.set_organism.text().strip()
        self._ruleset.description = self.set_description.text().strip()
        return self._ruleset

    def _problems(self) -> List[str]:
        return lib.validate(reg.to_dict(self._harvest()))

    def _save(self) -> None:
        self._write(self._path)

    def _save_as(self) -> None:
        self._write("")

    def _write(self, path: str) -> None:
        problems = self._problems()
        if problems:
            QMessageBox.warning(
                self, "Rule set not saved",
                "These problems must be fixed first:\n\n• " + "\n• ".join(problems[:12]))
            return
        ruleset = self._harvest()
        if not path:
            name, ok = QInputDialog.getText(self, "Save rule set", "Name:",
                                            text=ruleset.name)
            if not ok or not name.strip():
                return
            ruleset.name = name.strip()
            info = lib.store(ruleset, name=name.strip())
        else:
            info = lib.store(ruleset, path=path)
        self.saved_path = info.path
        QMessageBox.information(
            self, "Saved",
            f"'{info.name}' saved with {info.n_rules} rule(s).\n\n{info.path}")
        self.accept()
