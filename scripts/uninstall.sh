#!/usr/bin/env bash
#
# Remove what install.sh created: the virtual environment and the launcher.
#
#   ./scripts/uninstall.sh          remove the environment and launcher
#   ./scripts/uninstall.sh --all    also delete downloaded databases in ~/.gsm_toolbox
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOVE_DATA=0
[ "${1:-}" = "--all" ] && REMOVE_DATA=1

rm -rf "$HERE/.venv" && echo "removed .venv"
rm -rf "$HERE/GSM ToolBox.app" 2>/dev/null && echo "removed GSM ToolBox.app" || true
rm -f "$HOME/.local/share/applications/gsm-toolbox.desktop" 2>/dev/null \
    && echo "removed the desktop entry" || true

if [ "$REMOVE_DATA" = "1" ]; then
    # Kept by default: this folder holds databases that can take hours to rebuild, and
    # a reinstall reuses them. Removing it is opt-in for exactly that reason.
    rm -rf "$HOME/.gsm_toolbox" && echo "removed ~/.gsm_toolbox (cached databases)"
else
    echo "kept ~/.gsm_toolbox (cached databases) — pass --all to remove it too"
fi

echo "The source folder itself was not touched; delete it if you no longer want it."
