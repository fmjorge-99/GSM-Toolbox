# PyInstaller spec for GSM ToolBox.
# Build with:  pyinstaller gsm_toolbox.spec
#
# Produces a one-folder build under dist/GSM_ToolBox/. One-folder (rather than
# --onefile) starts much faster, which matters for the heavy scientific stack.

import os
import sys

from PyInstaller.utils.hooks import (
    collect_all,
    collect_data_files,
    collect_dynamic_libs,
    collect_submodules,
)

block_cipher = None


def _read_version(default="0.0.0"):
    """The package version, without importing the package (its deps may be absent)."""
    import re
    try:
        with open(os.path.join("gsm_toolbox", "__init__.py"), encoding="utf-8") as fh:
            match = re.search(r'^__version__\s*=\s*"([^"]+)"', fh.read(), re.M)
        return match.group(1) if match else default
    except OSError:
        return default


# --- Refuse to build while the previous build is still running --------
# PyInstaller cannot overwrite a running executable. It reports the PermissionError
# in passing and then *exits 0*, so the build looks successful while dist/ still holds
# the old binary — a stale exe that silently lacks every change just made. Fail here
# instead, with the instruction needed to fix it.
def _assert_exe_not_running():
    exe = os.path.join(os.path.abspath(SPECPATH), "dist", "GSM_ToolBox",
                       "GSM_ToolBox.exe")
    if not os.path.exists(exe):
        return
    try:
        # Opening for append needs the same write access the build will need.
        with open(exe, "ab"):
            pass
    except PermissionError:
        raise SystemExit(
            "Build aborted — GSM_ToolBox.exe is locked, which almost always means a "
            "copy of the app is still running.\n"
            "Close it (or: taskkill /IM GSM_ToolBox.exe /F) and build again.\n"
            "Building anyway would leave dist/ holding the PREVIOUS binary while "
            "appearing to succeed.")


_assert_exe_not_running()

# --- Import smoke-check (fail the build loudly on a broken module) ----
# A truncated/syntactically broken module would otherwise be packaged silently
# and only crash the user at runtime. Import every gsm_toolbox module now so
# `pyinstaller gsm_toolbox.spec` aborts before producing a non-functional binary.
def _smoke_import_all():
    import importlib
    import pkgutil
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    sys.path.insert(0, os.path.abspath(SPECPATH))
    import gsm_toolbox
    failures = []
    for mod in pkgutil.walk_packages(
            [os.path.dirname(gsm_toolbox.__file__)], prefix="gsm_toolbox."):
        if mod.name.endswith("__main__"):
            continue
        try:
            importlib.import_module(mod.name)
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{mod.name}: {exc}")
    if failures:
        raise SystemExit(
            "Build aborted — these modules failed to import:\n  "
            + "\n  ".join(failures))
    print("[spec] import smoke-check passed for all gsm_toolbox modules.")


_smoke_import_all()

# --- Anaconda DLL fix -------------------------------------------------
# This project's Python is Anaconda, which keeps the runtime DLLs that several
# stdlib C-extensions depend on in <base_prefix>/Library/bin rather than DLLs/.
# PyInstaller doesn't look there, so modules like _ctypes, _bz2, _lzma, _ssl,
# _hashlib and _sqlite3 fail at runtime with "DLL load failed". Bundle the whole
# family explicitly. Globs cover version-suffixed names (e.g. libcrypto-3-x64).
import glob as _glob

_extra_binaries = []
_lib_bin = os.path.join(sys.base_prefix, "Library", "bin")
if os.path.isdir(_lib_bin):
    _patterns = [
        "ffi*.dll",          # _ctypes
        "libbz2*.dll", "bz2*.dll",        # _bz2
        "liblzma*.dll", "lzma*.dll",      # _lzma
        "libssl*.dll", "libcrypto*.dll",  # _ssl, _hashlib
        "sqlite3*.dll",      # _sqlite3
        "zlib*.dll",         # zlib
        "liblz4*.dll", "libzstd*.dll",    # extras some builds link against
    ]
    _seen = set()
    for _pat in _patterns:
        for _p in _glob.glob(os.path.join(_lib_bin, _pat)):
            # Skip debug/trash artifacts; only ship real DLLs.
            if _p.lower().endswith(".dll") and _p not in _seen:
                _seen.add(_p)
                _extra_binaries.append((_p, "."))

# Bundle our example models + any future resources.
datas = [("gsm_toolbox/resources", "gsm_toolbox/resources")]
# cobra ships SBML schemas / example data it loads at runtime.
datas += collect_data_files("cobra")
# optlang/swiglpk solver backends are imported dynamically.
datas += collect_data_files("optlang")

datas += collect_data_files("straindesign")
# SCIP solver via pyscipopt (bundles the SCIP shared libraries).
datas += collect_data_files("pyscipopt")
_extra_binaries += collect_dynamic_libs("pyscipopt")

hiddenimports = []
# Package EVERY gsm_toolbox module explicitly. Many are imported lazily inside handler
# bodies (e.g. the preferences dialog, the enzyme/Selenzyme dialogs), which PyInstaller's
# static analysis can miss — and a missing GUI module silently drops the feature from the
# shipped app. This guarantees the whole application is bundled.
hiddenimports += collect_submodules("gsm_toolbox")
hiddenimports += collect_submodules("cobra")
hiddenimports += collect_submodules("optlang")
hiddenimports += collect_submodules("straindesign")
hiddenimports += collect_submodules("scipy")
# numpy 1.26 ships numpy/_core/ purely as a FORWARD-COMPATIBILITY shim: data pickled by
# numpy 2 refers to `numpy._core.multiarray`, and these modules redirect that to the
# numpy 1 `numpy.core`. Nothing imports them statically — pickle reaches them through
# find_class at load time — so PyInstaller's analysis never sees them and leaves them out.
# The result was an MDF analysis that worked in development and died in the packaged app
# with "No module named 'numpy._core'" while loading eQuilibrator's component-contribution
# parameters (a numpy 2-era .npz).
hiddenimports += collect_submodules("numpy")
hiddenimports += collect_submodules("pyscipopt")
hiddenimports += ["cdd"]  # pycddlib (EFM enumeration); imported lazily
_extra_binaries += collect_dynamic_libs("cdd")
hiddenimports += [
    "swiglpk", "scipy.optimize", "scipy.sparse.csgraph._validation",
    # matplotlib Qt backend used by the analysis plots
    "matplotlib.backends.backend_qtagg",
    # cloudpickle ships analysis closures to the worker process (workers.py)
    "cloudpickle",
    # QtWebEngine hosts the interactive Escher maps (Escher Explorer tab). Naming
    # these pulls in PyInstaller's PySide6 hooks, which bundle QtWebEngineProcess,
    # its resources/locales/ICU data and the QtWebChannel bridge (#T6).
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebChannel",
]
# openpyxl (+ et_xmlfile) lets the Omics dataset-prep tool read Excel tables (#F4);
# collect all submodules so pandas' lazy `import openpyxl.*` resolves in the frozen app.
hiddenimports += collect_submodules("openpyxl")
hiddenimports += ["et_xmlfile"]

# RDKit: clean 2-D structure drawings. collect_all grabs its compiled extensions,
# data files (atom/periodic-table data under RDConfig) and submodules.
_rdkit_datas, _rdkit_binaries, _rdkit_hidden = collect_all("rdkit")
datas += _rdkit_datas
_extra_binaries += _rdkit_binaries
hiddenimports += _rdkit_hidden

# eQuilibrator (thermodynamics / MDF, #2.1). The API + its cache/contribution packages
# and their data files (pint unit registry, group-decomposition tables) are pulled in via
# collect_all; the heavy compound cache is bundled separately below.
for _eqpkg in ("equilibrator_api", "equilibrator_cache", "component_contribution",
               "pint", "flexcache", "flexparser"):
    try:
        _d, _b, _h = collect_all(_eqpkg)
        datas += _d
        _extra_binaries += _b
        hiddenimports += _h
    except Exception as _e:  # noqa: BLE001 — thermo is optional; never fail the build
        print(f"[spec] eQuilibrator package '{_eqpkg}' not collected: {_e}")
hiddenimports += ["uncertainties", "sqlalchemy", "pooch", "slugify", "path"]

# The ~1.34 GB eQuilibrator compound cache is deliberately NOT bundled. The
# thermodynamics (MDF) suite is opt-in: it stays hidden until the user enables it in
# Settings ▸ Preferences ▸ "Enable MDF Suite", at which point the app asks for consent
# and downloads the dataset once. Shipping it by default tripled the installer size for a
# feature most users never touch, and presented thermodynamics as a headline capability
# when it is an advanced, assumption-laden analysis. Only the (small) eQuilibrator code
# packages above are bundled, so enabling the suite needs a data download and no reinstall.
print("[spec] eQuilibrator compound cache intentionally NOT bundled "
      "(MDF suite is opt-in; data downloads on user consent)")

a = Analysis(
    ["run_gsm_toolbox.py"],
    pathex=[],
    binaries=_extra_binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "PyQt5", "PyQt6"],
    cipher=block_cipher,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# --- Platform-specific packaging ------------------------------------------------
#
# PyInstaller does not cross-compile: each platform's bundle must be built on that
# platform. What differs between them is only the icon format and, on macOS, the extra
# .app wrapper — so the spec branches here rather than being forked three ways.
_ICONS = os.path.join("gsm_toolbox", "resources", "icons")
if sys.platform == "darwin":
    _icon = os.path.join(_ICONS, "app_icon.icns")
elif sys.platform == "win32":
    _icon = os.path.join(_ICONS, "app_icon.ico")
else:
    _icon = None                     # Linux takes its icon from the .desktop entry

if _icon and not os.path.exists(_icon):
    print(f"[spec] no icon at {_icon} — building without one")
    _icon = None

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GSM_ToolBox",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    # UPX mangles Mach-O binaries and macOS then refuses to load them; it also breaks
    # code signing. Off everywhere but Windows.
    upx=(sys.platform == "win32"),
    console=False,  # windowed app (no console)
    icon=_icon,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=(sys.platform == "win32"),
    name="GSM_ToolBox",
)

if sys.platform == "darwin":
    # Wraps the collected folder as a double-clickable .app. Without this, macOS users
    # get a bare Unix executable that Finder will not launch.
    app = BUNDLE(
        coll,
        name="GSM ToolBox.app",
        icon=_icon,
        bundle_identifier="io.github.gsmtoolbox",
        version=_read_version(),
        info_plist={
            "CFBundleName": "GSM ToolBox",
            "CFBundleDisplayName": "GSM ToolBox",
            "CFBundleShortVersionString": _read_version(),
            "NSHighResolutionCapable": True,
            # The app is a single-window desktop tool, not a document editor.
            "LSApplicationCategoryType": "public.app-category.education",
            # Qt reads this; without it macOS may start the app in the wrong appearance.
            "NSRequiresAquaSystemAppearance": False,
        },
    )
