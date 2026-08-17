#!/usr/bin/env bash
#
# Launch GSM ToolBox on Linux or macOS, installing on first run if needed.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="$HERE/.venv"

if [ ! -x "$VENV/bin/python" ]; then
    echo "No environment yet — running the installer first."
    "$HERE/scripts/install.sh" --no-desktop
fi

# A desktop application needs a display. Qt's own failure here is
# "could not connect to display", which sends people looking in the wrong place.
if [ "$(uname -s)" != "Darwin" ] \
   && [ -z "${DISPLAY:-}" ] && [ -z "${WAYLAND_DISPLAY:-}" ]; then
    echo "error: no DISPLAY or WAYLAND_DISPLAY — this is a desktop application." >&2
    echo "Over SSH use 'ssh -X'; headless, use: xvfb-run -a $0" >&2
    exit 1
fi

exec "$VENV/bin/python" "$HERE/run_gsm_toolbox.py" "$@"
