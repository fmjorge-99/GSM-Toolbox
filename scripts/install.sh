#!/usr/bin/env bash
#
# GSM ToolBox — install from source on Linux or macOS.
#
#   ./scripts/install.sh              install, and add a desktop/Launchpad entry
#   ./scripts/install.sh --no-desktop install only
#   ./scripts/install.sh --deps       print the system packages needed, then exit
#
# Creates a virtual environment in .venv beside this checkout, installs the
# dependencies into it, and registers a launcher so the application can be started
# without a terminal. Safe to re-run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="$HERE/.venv"
APP_NAME="GSM ToolBox"
WANT_DESKTOP=1

for arg in "$@"; do
    case "$arg" in
        --no-desktop) WANT_DESKTOP=0 ;;
        --deps) PRINT_DEPS=1 ;;
        -h|--help) sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $arg (try --help)" >&2; exit 2 ;;
    esac
done

case "$(uname -s)" in
    Darwin) OS=macos ;;
    Linux)  OS=linux ;;
    *) echo "error: this script covers Linux and macOS. On Windows use scripts\\run_windows.bat." >&2
       exit 1 ;;
esac

# ---------------------------------------------------------------------------
# System libraries. Qt fails at *runtime*, not at pip-install time, when these are
# missing — so the failure lands at first launch with an opaque message unless we
# get ahead of it.
# ---------------------------------------------------------------------------
linux_packages() {
    if command -v apt-get >/dev/null 2>&1; then
        local sound="libasound2t64"
        apt-cache show libasound2t64 >/dev/null 2>&1 || sound="libasound2"
        echo "sudo apt-get install -y python3-venv python3-pip libgl1 libegl1 \\"
        echo "  libxkbcommon-x11-0 libdbus-1-3 libnss3 libxdamage1 libxcomposite1 \\"
        echo "  libxrandr2 libxtst6 $sound libxcb-cursor0 libxcb-icccm4 \\"
        echo "  libxcb-keysyms1 libxcb-shape0 libxcb-xinerama0"
    elif command -v dnf >/dev/null 2>&1; then
        echo "sudo dnf install -y python3-virtualenv mesa-libGL mesa-libEGL \\"
        echo "  libxkbcommon-x11 dbus-libs nss libXdamage libXcomposite libXrandr \\"
        echo "  libXtst alsa-lib xcb-util-cursor"
    elif command -v pacman >/dev/null 2>&1; then
        echo "sudo pacman -S --needed python-virtualenv libgl libxkbcommon-x11 nss \\"
        echo "  libxdamage libxcomposite libxrandr libxtst alsa-lib xcb-util-cursor"
    elif command -v zypper >/dev/null 2>&1; then
        echo "sudo zypper install -y python3-virtualenv Mesa-libGL1 Mesa-libEGL1 \\"
        echo "  libxkbcommon-x11-0 libdbus-1-3 mozilla-nss libXdamage1 libXcomposite1 \\"
        echo "  libXrandr2 libXtst6 alsa xcb-util-cursor0"
    else
        echo "# Unrecognised distribution. Install the Qt/xcb runtime libraries listed"
        echo "# in INSTALL.md using your package manager."
    fi
}

if [ "${PRINT_DEPS:-0}" = "1" ]; then
    if [ "$OS" = macos ]; then
        echo "macOS needs no extra system libraries; Qt ships complete in the wheel."
        echo "You do need Python 3.10-3.12:   brew install \"python@3.12\""
    else
        linux_packages
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# Python. One dependency is pinned to the NumPy 1.x ABI, which has no 3.13 build.
# ---------------------------------------------------------------------------
find_python() {
    local candidate
    for candidate in python3.12 python3.11 python3.10 python3 python; do
        command -v "$candidate" >/dev/null 2>&1 || continue
        if "$candidate" -c 'import sys; raise SystemExit(
            0 if (3, 10) <= sys.version_info < (3, 13) else 1)' 2>/dev/null; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

if ! PYTHON="$(find_python)"; then
    echo "error: GSM ToolBox needs Python 3.10, 3.11 or 3.12." >&2
    echo "Found: $(python3 --version 2>&1 || echo 'no python3 on PATH')" >&2
    if [ "$OS" = macos ]; then
        echo 'Install one with:  brew install "python@3.12"' >&2
    else
        echo "Install one with your package manager, then re-run this script." >&2
    fi
    exit 1
fi
echo "==> Using $("$PYTHON" --version 2>&1) ($PYTHON)"

if [ ! -d "$VENV" ]; then
    echo "==> Creating virtual environment in .venv"
    "$PYTHON" -m venv "$VENV"
fi

echo "==> Installing dependencies (several minutes on a first run)"
"$VENV/bin/python" -m pip install --upgrade pip >/dev/null
"$VENV/bin/python" -m pip install -r "$HERE/requirements.txt"

# ---------------------------------------------------------------------------
# Verify before claiming success. A pip install that "worked" while Qt cannot
# initialise is the most common failure, and it is silent until first launch.
# ---------------------------------------------------------------------------
echo "==> Checking the installation"
if "$VENV/bin/python" - <<'PY'
import sys
try:
    import cobra                       # noqa: F401
    from PySide6.QtWidgets import QApplication   # noqa: F401
except Exception as exc:               # noqa: BLE001
    print(f"import failed: {exc}")
    sys.exit(1)
PY
then
    echo "    core libraries import cleanly"
else
    echo
    echo "warning: the dependencies installed but could not be imported." >&2
    if [ "$OS" = linux ]; then
        echo "You are probably missing system libraries. Install these and re-run:" >&2
        echo >&2
        linux_packages >&2
    fi
    exit 1
fi

# Qt can import and still fail to open a display; check that separately so the two
# problems are not confused with one another.
if ! QT_QPA_PLATFORM=offscreen "$VENV/bin/python" -c \
        "from PySide6.QtWidgets import QApplication; QApplication([])" >/dev/null 2>&1; then
    echo "    note: Qt could not start even offscreen — see INSTALL.md troubleshooting"
else
    echo "    Qt initialises"
fi

# ---------------------------------------------------------------------------
# Desktop integration
# ---------------------------------------------------------------------------
make_macos_app() {
    # A .app is a directory with a known layout; no macOS-only tooling is required
    # to create one. This wrapper runs the source install, so it stays correct when
    # the checkout is updated.
    local app="$HERE/GSM ToolBox.app"
    rm -rf "$app"
    mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources"

    cat > "$app/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
  "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key><string>${APP_NAME}</string>
    <key>CFBundleDisplayName</key><string>${APP_NAME}</string>
    <key>CFBundleIdentifier</key><string>io.github.gsmtoolbox</string>
    <key>CFBundleExecutable</key><string>gsm-toolbox</string>
    <key>CFBundleIconFile</key><string>app_icon.icns</string>
    <key>CFBundlePackageType</key><string>APPL</string>
    <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST

    cat > "$app/Contents/MacOS/gsm-toolbox" <<LAUNCH
#!/bin/bash
exec "$VENV/bin/python" "$HERE/run_gsm_toolbox.py" "\$@"
LAUNCH
    chmod +x "$app/Contents/MacOS/gsm-toolbox"

    local icns="$HERE/gsm_toolbox/resources/icons/app_icon.icns"
    [ -f "$icns" ] && cp "$icns" "$app/Contents/Resources/app_icon.icns"
    echo "    created $app — drag it to Applications if you like"
}

make_linux_entry() {
    local dir="$HOME/.local/share/applications"
    local icon="$HERE/gsm_toolbox/resources/icons/app_icon_256.png"
    mkdir -p "$dir"
    cat > "$dir/gsm-toolbox.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=${APP_NAME}
GenericName=Genome-scale metabolic modelling
Comment=Constraint-based metabolic modelling and strain design
Exec=$HERE/scripts/run.sh
Icon=$icon
Terminal=false
Categories=Science;Biology;Education;
Keywords=metabolism;FBA;COBRA;systems biology;
DESKTOP
    chmod +x "$dir/gsm-toolbox.desktop" 2>/dev/null || true
    command -v update-desktop-database >/dev/null 2>&1 \
        && update-desktop-database "$dir" >/dev/null 2>&1 || true
    echo "    added $dir/gsm-toolbox.desktop"
}

if [ "$WANT_DESKTOP" = "1" ]; then
    echo "==> Adding a launcher"
    if [ "$OS" = macos ]; then make_macos_app; else make_linux_entry; fi
fi

chmod +x "$HERE/scripts/"*.sh 2>/dev/null || true

cat <<DONE

Done.

  Launch:      ./scripts/run.sh
DONE
if [ "$WANT_DESKTOP" = "1" ]; then
    if [ "$OS" = macos ]; then
        echo "  or open:     \"$HERE/GSM ToolBox.app\""
    else
        echo "  or find \"$APP_NAME\" in your applications menu"
    fi
fi
echo "  Uninstall:   ./scripts/uninstall.sh"
echo
