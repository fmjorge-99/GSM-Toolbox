"""GSM ToolBox application entry point."""

from __future__ import annotations

import sys


def main() -> int:
    # A windowed build has no console, so sys.stdout/sys.stderr are None. Libraries that
    # write to them regardless — tqdm, which eQuilibrator uses — crash the app with
    # "'NoneType' object has no attribute 'write'". Fix this FIRST: it must be in place
    # before any dependency is imported, let alone run.
    from .core import numpy_compat, stdio_guard

    stdio_guard.install()
    # Data pickled by numpy 2 (eQuilibrator's parameter file) must remain loadable under
    # the numpy 1 this app pins. Aliases have to exist before anything unpickles.
    numpy_compat.install()

    # Guard against multiprocessing workers (e.g. from StrainDesign) re-launching
    # the GUI when frozen — without this, spawned workers open extra windows.
    import multiprocessing

    multiprocessing.freeze_support()

    # Import Qt lazily so that ``import gsm_toolbox`` stays GUI-free.
    import os

    from PySide6.QtGui import QIcon
    from PySide6.QtWidgets import QApplication

    from . import __app_name__
    from .gui.main_window import MainWindow
    from .gui.style import STYLESHEET
    from .resources import app_icon_path

    app = QApplication(sys.argv)
    app.setApplicationName(__app_name__)
    app.setStyleSheet(STYLESHEET)
    icon_file = app_icon_path()
    if os.path.exists(icon_file):
        app.setWindowIcon(QIcon(icon_file))
    window = MainWindow()
    window.showMaximized()

    # Optional self-test: trigger a StrainDesign run (which spawns a multiprocessing
    # pool) to verify the frozen build does not open extra windows. Enabled only via
    # the GSM_SELFTEST env var; a no-op in normal use.
    import os

    if os.environ.get("GSM_SELFTEST") == "thermo":
        from PySide6.QtCore import QTimer

        def _selftest():
            # The packaged-app failure this guards against: eQuilibrator's parameter file
            # is a numpy 2-era pickle, and unpickling it needs `numpy._core` — a shim no
            # code imports statically, so a bundler can silently omit it. Proving the
            # import resolves here is proving the MDF analysis can start at all.
            import io
            import pickle
            ok = True
            try:
                pickle.Unpickler(io.BytesIO(b".")).find_class(
                    "numpy._core.multiarray", "_reconstruct")
                print("NUMPY_CORE_SELFTEST resolved=True", flush=True)
            except Exception as exc:  # noqa: BLE001
                ok = False
                print(f"NUMPY_CORE_SELFTEST resolved=False error={exc}", flush=True)
            try:
                from .core import thermodynamics as thermo
                n = len(thermo._accessions(_probe_metabolite()))
                print(f"THERMO_SELFTEST accessions={n}", flush=True)
                ok = ok and n >= 3
            except Exception as exc:  # noqa: BLE001
                ok = False
                print(f"THERMO_SELFTEST failed error={exc}", flush=True)
            app.exit(0 if ok else 5)

        QTimer.singleShot(1500, _selftest)

    if os.environ.get("GSM_SELFTEST") == "rdkit":
        from PySide6.QtCore import QTimer

        def _selftest():
            # Verify RDKit (bundled) can render a clean 2-D structure in the frozen build.
            from .gui.widgets.structure_fetcher import _render_bw
            png = _render_bw(smiles="CC(=O)C(=O)[O-]", size=200)
            ok = bool(png) and png[:4] == b"\x89PNG"
            print(f"RDKIT_SELFTEST bytes={len(png)} ok={ok}", flush=True)
            app.exit(0 if ok else 5)

        QTimer.singleShot(1200, _selftest)

    if os.environ.get("GSM_SELFTEST") == "straindesign":
        from PySide6.QtCore import QTimer

        def _selftest():
            from .core import io_models
            from .core.analysis import strain_design as sd
            from .resources import example_model_path

            model = io_models.load_model(example_model_path())
            model.reactions.EX_o2_e.lower_bound = 0
            try:
                sd.compute_knockout_design(model, "EX_etoh_e", max_knockouts=2,
                                           max_solutions=2, time_limit=40)
            except Exception:  # noqa: BLE001
                pass
            app.quit()

        QTimer.singleShot(1500, _selftest)

    # Optional self-test: the restructured interface, checked *in the frozen build*.
    #
    # Worth running here rather than only under pytest because the failure it guards
    # against is a lifetime bug, not a logic one: Qt menus created by ``addMenu(title)``
    # are owned by Python, and a packaged build collects on a different schedule than a
    # test process. A menu bar that lists menus whose C++ objects have been freed looks
    # perfectly normal until something touches one.
    if os.environ.get("GSM_SELFTEST") == "ui":
        import gc

        from PySide6.QtCore import QTimer

        def _selftest():
            problems = []
            gc.collect()
            for menu in window._menus:
                try:
                    menu.actions()
                except RuntimeError:
                    problems.append(f"menu freed: {menu!r}")
            # Reaching a menu through the bar makes a second wrapper — the path that
            # actually triggered the bug.
            for action in window.menuBar().actions():
                sub = action.menu()
                if sub is None:
                    continue
                try:
                    sub.actions()
                except RuntimeError:
                    problems.append(f"menu freed via menuBar(): {action.text()}")
            gc.collect()
            try:
                names = [m.title() for m in window._menus]
            except RuntimeError:
                names = []
                problems.append("menus freed after a second collection")

            tabs = [window.tabs.tabText(i) for i in range(window.tabs.count())]
            if tabs != ["Analysis", "Pathway Design", "Dynamic Analysis", "Network Map"]:
                problems.append(f"unexpected tabs: {tabs}")
            if getattr(window, "_network_button", None) is None:
                problems.append("no Network Visualization button on the toolbar")
            if window.act_media.isEnabled():
                problems.append("model-dependent action enabled with no model loaded")

            # Rule sets ship as importable examples, never as an applied default. If the
            # resource were missing from the bundle the library would look empty and the
            # Synechocystis set would be unreachable from the packaged app.
            from .core import rule_library as rules
            examples = rules.examples()
            if not examples:
                problems.append("no example rule sets bundled")
            for _label, path in examples:
                info = rules.inspect(path)
                if not info.valid:
                    problems.append(f"bundled example invalid: {info.problems[:2]}")

            panel = getattr(window, "dynamic_panel", None)
            for attribute in ("rules_btn", "use_regulation", "substrates",
                              "scan_targets_btn", "plot_btn"):
                if getattr(panel, attribute, None) is None:
                    problems.append(f"dynamic panel missing {attribute}")

            print(f"UI_SELFTEST menus={len(names)} tabs={len(tabs)} "
                  f"examples={len(examples)} problems={problems or 'none'}", flush=True)
            app.exit(0 if not problems else 6)

        QTimer.singleShot(1500, _selftest)

    # Optional self-test: exercise the Omics dataset-prep pipeline including an Excel
    # round-trip — verifies openpyxl is bundled in the frozen build (#F4).
    if os.environ.get("GSM_SELFTEST") == "omics":
        import tempfile

        import pandas as pd

        from .core import io_models, omics_prep
        from .resources import example_model_path
        model = io_models.load_model(example_model_path())
        genes = [g.id for g in model.genes][:5]
        xlsx = os.path.join(tempfile.gettempdir(), "gsm_omics_selftest.xlsx")
        pd.DataFrame({"gene": genes, "s1": range(5), "s2": range(5, 10)}).to_excel(xlsx, index=False)
        df = omics_prep.read_table(xlsx)   # needs openpyxl to read .xlsx
        ds = omics_prep.prepare_dataset(df, model, kind="transcriptomics",
                                        id_column="gene", value_columns=["s1", "s2"])
        ok = ds.summary.n_model_targets == 5
        print(f"OMICS_SELFTEST mapped={ds.summary.n_model_targets} ok={ok}", flush=True)
        return 0 if ok else 7

    # Optional self-test: render an interactive Escher map in the embedded
    # QtWebEngine view — verifies the frozen build bundled QtWebEngineProcess, its
    # resources and the escher assets, and that a map draws (circles > 0) (#T6).
    if os.environ.get("GSM_SELFTEST") == "escher":
        from PySide6.QtCore import QTimer

        from .core import io_models
        from .resources import example_model_path

        def _selftest():
            model = io_models.load_model(example_model_path())
            panel = window.escher_explorer
            panel.set_model(model)
            panel.rebuild()

            def _probe():
                def _done(val):
                    import json as _json
                    try:
                        d = _json.loads(val) if val else {}
                    except Exception:  # noqa: BLE001
                        d = {}
                    circles = int(d.get("circles", 0))
                    ok = circles > 0
                    print(f"ESCHER_SELFTEST circles={circles} ok={ok}", flush=True)
                    app.exit(0 if ok else 6)
                panel.view.page().runJavaScript(
                    "JSON.stringify({circles:document.querySelectorAll('#map circle').length})",
                    0, _done)

            QTimer.singleShot(4000, _probe)

        QTimer.singleShot(2000, _selftest)

    # Optional self-test: run a job through the out-of-process worker (workers.py)
    # to verify the frozen build can spawn/return from the analysis subprocess.
    if os.environ.get("GSM_SELFTEST") == "worker":
        from PySide6.QtCore import QTimer

        from .gui.workers import JobRunner

        def _worker_job():
            import os as _os

            from .core import io_models
            from .core.analysis import fba
            from .resources import example_model_path
            model = io_models.load_model(example_model_path())
            return _os.getpid(), fba.run_fba(model).objective_value

        def _selftest():
            runner = JobRunner(window)

            def _done(result):
                pid, val = result
                ok = pid != os.getpid() and abs(val - 0.8739) < 1e-2
                print(f"WORKER_SELFTEST pid={pid} parent={os.getpid()} "
                      f"obj={val:.4f} out_of_process={pid != os.getpid()} ok={ok}",
                      flush=True)
                app.exit(0 if ok else 3)

            runner.done.connect(_done)
            runner.error.connect(lambda e: (print("WORKER_SELFTEST error:", e, flush=True),
                                            app.exit(4)))
            window._selftest_runner = runner  # keep alive
            runner.start(_worker_job)

        QTimer.singleShot(1500, _selftest)

    return app.exec()


def _probe_metabolite():
    """A metabolite annotated the way the bundled database writes them.

    Used by the ``thermo`` self-test: the identifier reader must cope with
    ``"KEGG Compound"``-style keys holding identifiers.org URLs, not just MIRIAM
    spellings, or the thermodynamics analysis silently resolves nothing.
    """
    import cobra

    met = cobra.Metabolite("btal_c", formula="C4H8O", compartment="c")
    met.name = "Butanal C4H8O"
    met.annotation = {
        "CHEBI": ["http://identifiers.org/chebi/CHEBI:13923"],
        "InChI Key": "https://identifiers.org/inchikey/ZTQSAGDEMFDKMZ-UHFFFAOYSA-N",
        "KEGG Compound": "http://identifiers.org/kegg.compound/C01412",
    }
    return met


if __name__ == "__main__":
    sys.exit(main())
