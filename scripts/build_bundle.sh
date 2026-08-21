#!/usr/bin/env bash
#
# Build a self-contained GSM ToolBox bundle for the machine you run this on.
#
#   macOS  ->  dist/GSM ToolBox.app  and  GSM_ToolBox-<version>-macOS-<arch>.dmg
#   Linux  ->  dist/GSM_ToolBox/     and  GSM_ToolBox-<version>-linux-<arch>.tar.gz
#
# PyInstaller does NOT cross-compile. A macOS bundle must be built on macOS and a
# Linux bundle on Linux; there is no way around this, which is why the release
# workflow in .github/workflows builds each one on its own runner.
#
# The result contains Python, Qt, the solvers and every dependency, so the user needs
# nothing installed. It is large (roughly 400-900 MB unpacked) for that reason.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"

case "$(uname -s)" in
    Darwin) OS=macos ;;
    Linux)  OS=linux ;;
    *) echo "error: run this on macOS or Linux. On Windows use the Inno Setup script." >&2
       exit 1 ;;
esac
ARCH="$(uname -m)"
VERSION="$(sed -n 's/^__version__ *= *"\([^"]*\)".*/\1/p' gsm_toolbox/__init__.py)"
[ -n "$VERSION" ] || { echo "error: could not read the version" >&2; exit 1; }
echo "==> GSM ToolBox $VERSION for $OS/$ARCH"

# Build in an isolated environment so a stale or partial dev venv cannot leak into
# the bundle — PyInstaller freezes whatever it finds importable.
BUILD_VENV="$HERE/.venv-build"
if [ ! -x "$BUILD_VENV/bin/python" ]; then
    for candidate in python3.12 python3.11 python3.10 python3; do
        command -v "$candidate" >/dev/null 2>&1 || continue
        "$candidate" -c 'import sys; raise SystemExit(
            0 if (3, 10) <= sys.version_info < (3, 13) else 1)' 2>/dev/null || continue
        echo "==> Creating build environment with $candidate"
        "$candidate" -m venv "$BUILD_VENV"
        break
    done
fi
[ -x "$BUILD_VENV/bin/python" ] || {
    echo "error: need Python 3.10-3.12 to build." >&2; exit 1; }

PY="$BUILD_VENV/bin/python"
"$PY" -m pip install --upgrade pip >/dev/null
echo "==> Installing dependencies and PyInstaller"
"$PY" -m pip install -r requirements.txt pyinstaller >/dev/null

echo "==> Freezing"
rm -rf build dist
"$PY" -m PyInstaller gsm_toolbox.spec --noconfirm --clean

# ---------------------------------------------------------------------------
if [ "$OS" = macos ]; then
    APP="dist/GSM ToolBox.app"
    [ -d "$APP" ] || { echo "error: the spec did not produce $APP" >&2; exit 1; }

    # Unsigned bundles are quarantined by Gatekeeper on download. An ad-hoc signature
    # does not avoid that, but it does stop the "damaged and can't be opened" error
    # that an unsigned, modified bundle otherwise triggers on Apple Silicon.
    if command -v codesign >/dev/null 2>&1; then
        echo "==> Ad-hoc signing"
        codesign --force --deep --sign - "$APP" 2>/dev/null \
            && echo "    signed" || echo "    signing skipped (not fatal)"
    fi

    DMG="GSM_ToolBox-${VERSION}-macOS-${ARCH}.dmg"
    echo "==> Building $DMG"
    rm -f "$DMG"
    STAGE="$(mktemp -d)"
    cp -R "$APP" "$STAGE/"
    ln -s /Applications "$STAGE/Applications"      # drag-to-install target
    hdiutil create -volname "GSM ToolBox $VERSION" -srcfolder "$STAGE" \
        -ov -format UDZO "$DMG" >/dev/null
    rm -rf "$STAGE"
    echo
    echo "Built:  $APP"
    echo "        $DMG  ($(du -h "$DMG" | cut -f1))"
    echo
    echo "Note: the bundle is unsigned and not notarised, so on first launch macOS"
    echo "will refuse it. Right-click the app and choose Open, or run:"
    echo "  xattr -dr com.apple.quarantine \"/Applications/GSM ToolBox.app\""
else
    OUT="dist/GSM_ToolBox"
    [ -d "$OUT" ] || { echo "error: the spec did not produce $OUT" >&2; exit 1; }

    # A .desktop file inside the tarball, with a relative Exec so it works wherever
    # the user unpacks it once they run install-desktop.sh.
    cat > "$OUT/gsm-toolbox.desktop" <<'DESKTOP'
[Desktop Entry]
Type=Application
Name=GSM ToolBox
GenericName=Genome-scale metabolic modelling
Comment=Constraint-based metabolic modelling and strain design
Exec=GSM_ToolBox
Icon=gsm-toolbox
Terminal=false
Categories=Science;Biology;Education;
DESKTOP
    cp gsm_toolbox/resources/icons/app_icon_256.png "$OUT/gsm-toolbox.png" 2>/dev/null || true

    cat > "$OUT/install-desktop.sh" <<'INSTALLER'
#!/usr/bin/env bash
# Register this unpacked bundle with the desktop environment.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APPS="$HOME/.local/share/applications"
ICONS="$HOME/.local/share/icons/hicolor/256x256/apps"
mkdir -p "$APPS" "$ICONS"
sed "s|^Exec=.*|Exec=$HERE/GSM_ToolBox|" "$HERE/gsm-toolbox.desktop" \
    > "$APPS/gsm-toolbox.desktop"
cp "$HERE/gsm-toolbox.png" "$ICONS/gsm-toolbox.png" 2>/dev/null || true
command -v update-desktop-database >/dev/null 2>&1 \
    && update-desktop-database "$APPS" >/dev/null 2>&1 || true
echo "GSM ToolBox added to your applications menu."
INSTALLER
    chmod +x "$OUT/install-desktop.sh"

    TARBALL="GSM_ToolBox-${VERSION}-linux-${ARCH}.tar.gz"
    echo "==> Building $TARBALL"
    rm -f "$TARBALL"
    tar -czf "$TARBALL" -C dist GSM_ToolBox
    echo
    echo "Built:  $OUT/"
    echo "        $TARBALL  ($(du -h "$TARBALL" | cut -f1))"
    echo
    echo "To use: tar -xzf $TARBALL && ./GSM_ToolBox/GSM_ToolBox"
    echo "        (optionally ./GSM_ToolBox/install-desktop.sh to add a menu entry)"
fi
