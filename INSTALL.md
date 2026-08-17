# Installation guide

| Route | Platform | Needs Python? |
|---|---|---|
| [Ready-made bundle](#1-ready-made-bundles) | Windows, macOS, Linux | **No** |
| [One-command install from source](#2-one-command-install-linux-and-macos) | Linux, macOS | Yes, 3.10–3.12 |
| [Manual install from source](#3-manual-install-any-platform) | Any | Yes, 3.10–3.12 |
| [Build a bundle yourself](#4-building-a-bundle-yourself) | Any | Yes |

**In a hurry?** On Linux or macOS:

```sh
./scripts/install.sh && ./scripts/run.sh
```

---

## 1. Ready-made bundles

Download from the [Releases](../../releases) page. Each bundle contains Python, Qt,
the solvers and every dependency, so nothing else needs installing.

### Windows, installer

Run `GSM_ToolBox_Setup_<version>.exe`. It installs per user into
`%LOCALAPPDATA%\Programs\GSM ToolBox`, so no administrator rights are needed. Launch it
from the Start menu afterwards.

Uninstall from Settings, Apps, or with the Start menu shortcut.

### Windows, portable

Download `GSM_ToolBox-<version>-windows-portable.zip`, unzip it anywhere, and run
`GSM_ToolBox.exe` inside the folder. Nothing is installed and nothing is written to the
registry, so this works from a USB stick or a shared drive, and on a machine where you
cannot install software.

Downloaded databases still go to your user profile at `%USERPROFILE%\.gsm_toolbox`. To
keep those with the portable copy instead, create an empty folder named `.gsm_toolbox`
beside `GSM_ToolBox.exe` and the application will use that.

> **SmartScreen.** Both builds are unsigned, so Windows may show "Windows protected your
> PC" the first time. Choose More info, then Run anyway. Signing requires a paid
> certificate.

### macOS

Open the `.dmg` and drag **GSM ToolBox** to Applications. Pick the file matching your
machine. Use `arm64` for Apple Silicon and `x86_64` for Intel.

> **First launch.** The app is unsigned and not notarised, so macOS will refuse it with
> *"cannot be opened because the developer cannot be verified"*. Either **right-click the
> app and choose Open** (then confirm), or clear the quarantine flag once:
> ```sh
> xattr -dr com.apple.quarantine "/Applications/GSM ToolBox.app"
> ```
> This happens because the build is not signed with a paid Apple Developer
> certificate. It is not caused by anything the application does.

### Linux

```sh
tar -xzf GSM_ToolBox-<version>-linux-x86_64.tar.gz
./GSM_ToolBox/GSM_ToolBox
```

Optionally add it to your applications menu:

```sh
./GSM_ToolBox/install-desktop.sh
```

The bundle is built on Ubuntu 22.04, so it needs glibc 2.35 or newer. That covers
Ubuntu 22.04 and later, Debian 12 and later, and Fedora 36 and later. On an older
distribution, install from source instead.

---

## 2. One-command install (Linux and macOS)

Installs into a virtual environment beside the checkout and registers a launcher, so the
application can be started without a terminal afterwards.

```sh
git clone <repository-url>
cd GSM_Toolbox_Distribution
./scripts/install.sh
```

The script picks a suitable Python, creates `.venv`, installs the dependencies,
**verifies that the libraries actually import and that Qt can initialise**, and then adds
a launcher:

- On Linux it adds an entry to your applications menu.
- On macOS it creates a double-clickable `GSM ToolBox.app` in the folder, which you can
  drag to Applications.

| Command | What it does |
|---|---|
| `./scripts/install.sh` | Install and add a launcher |
| `./scripts/install.sh --no-desktop` | Install only |
| `./scripts/install.sh --deps` | Print the system packages you need, then exit |
| `./scripts/run.sh` | Launch (installs first if needed) |
| `./scripts/uninstall.sh` | Remove the environment and launcher |
| `./scripts/uninstall.sh --all` | Also delete cached databases in `~/.gsm_toolbox` |

### Linux system libraries

Qt needs shared libraries that are missing from minimal installs. `install.sh` detects
your package manager and prints the exact command. Run `./scripts/install.sh --deps` to
see it without installing anything. For reference:

**Debian / Ubuntu**
```sh
sudo apt install -y python3-venv python3-pip libgl1 libegl1 libxkbcommon-x11-0 \
  libdbus-1-3 libnss3 libxdamage1 libxcomposite1 libxrandr2 libxtst6 libasound2t64 \
  libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1 libxcb-shape0 libxcb-xinerama0
```
On releases older than Ubuntu 24.04, use `libasound2` instead of `libasound2t64`.

**Fedora / RHEL**
```sh
sudo dnf install -y python3-virtualenv mesa-libGL mesa-libEGL libxkbcommon-x11 \
  dbus-libs nss libXdamage libXcomposite libXrandr libXtst alsa-lib xcb-util-cursor
```

**Arch**
```sh
sudo pacman -S --needed python-virtualenv libgl libxkbcommon-x11 nss \
  libxdamage libxcomposite libxrandr libxtst alsa-lib xcb-util-cursor
```

macOS needs no extra system libraries, because Qt ships complete in the wheel. You do
need Python 3.10 to 3.12: `brew install "python@3.12"`.

---

## 3. Manual install (any platform)

If you would rather not use the helper script, or you are on Windows and want to
run from source:

### Install and run

```sh
git clone <repository-url>
cd GSM_Toolbox_Distribution

python3 -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate

pip install --upgrade pip
pip install -r requirements.txt

python run_gsm_toolbox.py
```

On Windows, `scripts\run_windows.bat` creates the environment and launches in one step.

### Headless servers and WSL

The application is a desktop GUI and needs a display.

- WSL2 on Windows 11 runs it directly through WSLg, with no extra setup.
- **A remote server** needs X forwarding (`ssh -X`) or a virtual display:
  ```sh
  sudo apt install -y xvfb
  xvfb-run -a python run_gsm_toolbox.py
  ```
- To confirm Qt works before anything else:
  ```sh
  python -c "from PySide6.QtWidgets import QApplication; QApplication([]); print('Qt OK')"
  ```

### macOS

```sh
brew install "python@3.12"
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python run_gsm_toolbox.py
```

Apple Silicon is supported. On first launch macOS may ask permission for network
access. The application uses it only to download reaction databases when you ask it to.

---

## 4. Building a bundle yourself

PyInstaller does not cross-compile. Each bundle must be built on the platform it is
for. There is no way to produce a macOS app from Windows, or the reverse.

### macOS and Linux

```sh
./scripts/build_bundle.sh
```

Builds in an isolated environment and produces a `.dmg` on macOS or a `.tar.gz` on
Linux, both self-contained.

### Windows

Requires [Inno Setup 6](https://jrsoftware.org/isinfo.php).

```powershell
pip install -r requirements.txt pyinstaller
python -m PyInstaller gsm_toolbox.spec --noconfirm
& "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe" installer\gsm_toolbox.iss
```

The frozen application lands in `dist\GSM_ToolBox\`; the installer in
`installer\Output\`. See [`installer/README.md`](installer/README.md).

### All three at once, without owning the machines

`.github/workflows/release.yml` builds Windows, both macOS architectures and Linux on
GitHub's runners and attaches the results to a draft release:

```sh
git tag v0.3.8 && git push origin v0.3.8
```

This is the practical route if you do not have a Mac. It is free for public
repositories.

---

## Where data is stored

Everything the application downloads or remembers lives under one folder:

| Platform | Location |
|---|---|
| Windows | `%USERPROFILE%\.gsm_toolbox` |
| Linux / macOS | `~/.gsm_toolbox` |

It holds cached reaction databases, molecule structure images, regulatory rule sets and
preferences. Deleting it resets the application; anything needed is fetched again.

Nothing is transmitted anywhere. Network access is outbound only, to the database
providers, and only when you ask for a download.

---

## Troubleshooting

**`could not load the Qt platform plugin "xcb"` (Linux)**
Install the packages listed above. To find which one is missing:
```sh
QT_DEBUG_PLUGINS=1 python run_gsm_toolbox.py 2>&1 | grep -i "cannot load\|undefined symbol"
```

**`ImportError: numpy.core.multiarray failed to import`**
NumPy 2.x was installed. Run `pip install "numpy>=1.24,<2"`.

**Python 3.13 install fails**
Use 3.10–3.12; a dependency is pinned to the NumPy 1.x ABI.

**Solver errors on strain design**
`pyscipopt` failed to install. The application falls back to GLPK, which is slower on
mixed-integer problems but correct. To retry: `pip install pyscipopt`.

**The application starts but analyses fail immediately**
Confirm the solver stack independently:
```sh
python -c "import cobra; m = cobra.io.load_model('textbook'); print(m.optimize())"
```
A number here means COBRApy and its solver are working and the problem lies elsewhere.
Please report it with that output.

**Windows: "DLL load failed while importing QtCore"**
The Microsoft Visual C++ Redistributable is missing or outdated. Install the current
x64 redistributable from Microsoft.
