# Installation guide

| Route | Platform | Needs Python? |
|---|---|---|
| [Download a ready-made build](#1-ready-made-builds) | **Windows only, for now** | No |
| [One command from source](#2-one-command-install-linux-and-macos) | macOS, Linux | Yes, 3.10 to 3.12 |
| [Manual install from source](#3-manual-install-any-platform) | Any | Yes, 3.10 to 3.12 |
| [Build a bundle yourself](#4-building-a-bundle-yourself) | Any | Yes |

> **macOS users, read this first.** The current release has Windows files only, so there
> is no `.dmg` to download yet. Installing from source takes about five minutes and is
> fully supported — jump to [macOS](#macos). Do not follow the Windows sections.

The quickest way to run the toolbox on any platform is from source:

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_gsm_toolbox.py
```

---

## 1. Ready-made builds

A frozen build contains Python, Qt and the solvers, so nothing else needs installing.
These files are around 230 to 350 MB, which is far above the size a git repository
accepts, so they are **not in the repository**. They are published on the
[Releases](../../releases) page when a release is made.

If the Releases page is empty, or has no file for your platform, build one yourself.
Section 4 covers Windows and takes about ten minutes. `scripts/build_bundle.sh` does the
same on macOS and Linux.

### Windows, installer

Run `GSM_ToolBox_Setup_<version>.exe`. It installs per user into
`%LOCALAPPDATA%\Programs\GSM ToolBox`, so no administrator rights are needed. Launch it
from the Start menu afterwards, and uninstall from Settings, Apps.

### Windows, portable

Unzip `GSM_ToolBox-<version>-windows-portable.zip` anywhere and run `GSM_ToolBox.exe`
inside the folder. Nothing is installed and nothing is written to the registry, so this
works from a USB stick and on a machine where you cannot install software.

Downloaded databases go to `%USERPROFILE%\.gsm_toolbox`. To keep them inside the portable
folder instead, create an empty folder named `.gsm_toolbox` next to `GSM_ToolBox.exe`
before starting the application.

### macOS

**There is no macOS download yet.** The published release contains Windows files only.
Install from source instead — see [macOS](#macos) under section 2, which is a single
command and takes about five minutes.

Once a `.dmg` is published, or once you have built one yourself with
`scripts/build_bundle.sh`, install it by opening the `.dmg` and dragging the app to
Applications. Use the `arm64` file on Apple Silicon and the `x86_64` file on Intel.

> **First launch of a downloaded build.** It is not signed with a paid Apple
> certificate, so macOS refuses it with "cannot be opened because the developer cannot be
> verified". Right-click the app and choose Open, then confirm. Alternatively clear the
> quarantine flag once:
> ```sh
> xattr -dr com.apple.quarantine "/Applications/GSM ToolBox.app"
> ```
> This applies only to a bundle you downloaded. A bundle you built yourself, and the
> source install below, are never quarantined and need none of this.

### Linux

```sh
tar -xzf GSM_ToolBox-<version>-linux-x86_64.tar.gz
./GSM_ToolBox/GSM_ToolBox
```

Add it to your applications menu with `./GSM_ToolBox/install-desktop.sh`.

The bundle is built on Ubuntu 22.04, so it needs glibc 2.35 or newer. That covers Ubuntu
22.04 and later, Debian 12 and later, and Fedora 36 and later. On an older distribution,
install from source instead.

> **SmartScreen.** The Windows builds are unsigned, so Windows may show "Windows
> protected your PC" the first time. Choose More info, then Run anyway. Signing requires
> a paid certificate.

---

## 2. One-command install (macOS and Linux)

Installs into a virtual environment beside the checkout and registers a launcher, so the
application can be started without a terminal afterwards.

> **Run the scripts with `bash scripts/...`, not `./scripts/...`.** The repository stores
> them without the executable bit, so `./scripts/install.sh` fails with `permission
> denied`. Prefixing with `bash` works regardless. If you prefer the short form, run
> `chmod +x scripts/*.sh` once first.

### macOS

**Step 1 — install Python 3.10, 3.11 or 3.12.** macOS ships Python 3.9, which is too old,
so this step is genuinely required. Either option is fine:

- **Without a package manager** (simplest): download the macOS 64-bit universal2
  installer for Python 3.12 from [python.org/downloads/macos](https://www.python.org/downloads/macos/)
  and double-click it. Nothing else to configure.
- **With Homebrew**: `brew install "python@3.12"`.

**Step 2 — install the toolbox.** Open Terminal and run:

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
bash scripts/install.sh
```

It takes a few minutes. The script finds your Python, creates `.venv`, installs the
dependencies, **checks that they really import and that Qt can start**, and then creates a
double-clickable **`GSM ToolBox.app`** in the folder.

**Step 3 — launch it.** Drag `GSM ToolBox.app` to your Applications folder and open it
like any other app, or run `bash scripts/run.sh`.

Apple Silicon and Intel are both supported, and no extra system libraries are needed
because Qt ships complete in the wheel. Gatekeeper does not interfere with this route:
nothing was downloaded as a packaged app, so nothing is quarantined. On first launch
macOS may ask for permission for network access, which the application uses only to
download reaction databases when you ask it to.

**If something goes wrong**

| Symptom | Cause and fix |
|---|---|
| `permission denied` | Use `bash scripts/install.sh`, or `chmod +x scripts/*.sh` first |
| `error: GSM ToolBox needs Python 3.10, 3.11 or 3.12` | Step 1 was skipped, or a newer Python is first on `PATH`. Install 3.12 and re-run |
| `xcrun: error: invalid active developer path` | `git` needs Apple's command line tools: run `xcode-select --install`, then retry |
| Build errors while pip installs | Almost always Python 3.13. Check with `python3 --version` and install 3.12 |

### Linux

```sh
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox
bash scripts/install.sh
```

The script detects your package manager and prints the exact system packages Qt needs if
any are missing, then adds an entry to your applications menu.

### All the commands

| Command | What it does |
|---|---|
| `bash scripts/install.sh` | Install and add a launcher |
| `bash scripts/install.sh --no-desktop` | Install only |
| `bash scripts/install.sh --deps` | Print the system packages you need, then exit |
| `bash scripts/run.sh` | Launch (installs first if needed) |
| `bash scripts/uninstall.sh` | Remove the environment and launcher |
| `bash scripts/uninstall.sh --all` | Also delete cached databases in `~/.gsm_toolbox` |

### Linux system libraries

Qt needs shared libraries that are missing from minimal installs. `install.sh` detects
your package manager and prints the exact command. Run `bash scripts/install.sh --deps` to
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
git clone https://github.com/fmjorge-99/GSM-Toolbox.git
cd GSM-Toolbox

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
bash scripts/build_bundle.sh
```

Builds in an isolated environment and produces a `.dmg` on macOS or a `.tar.gz` on
Linux, both self-contained.

**On macOS this is the second way to install**, and the one to use if you want a proper
app in Applications that does not depend on the checkout, or if you need to hand a
single file to colleagues who have no Python. It needs Python 3.10 to 3.12 once, to
build with, and takes roughly fifteen minutes. It produces:

- `dist/GSM ToolBox.app` — drag it to Applications and it is done
- `GSM_ToolBox-<version>-macOS-<arch>.dmg` — the file to share

The script ad-hoc signs the app, which is what prevents the "damaged and can't be
opened" error Apple Silicon otherwise shows for an unsigned bundle. The `.dmg` is still
unsigned and unnotarised, so anyone who **downloads** it will need the right-click ▸ Open
step described in section 1. The app you built locally is not quarantined and opens
normally.

### Windows

You need [Inno Setup 6](https://jrsoftware.org/isinfo.php) and Python 3.10 to 3.12.

This is two steps in a fixed order. PyInstaller freezes the application into
`dist\GSM_ToolBox`, and Inno Setup then wraps that folder into a single setup file. The
Inno script only packages what already exists, so running it first fails with
"No files found matching ...\dist\GSM_ToolBox\*". That message means the freeze step has
not run, not that a file is missing from the repository.

```powershell
# from the repository root
pip install -r requirements.txt pyinstaller

# step 1, about 8 minutes: freeze the application
python -m PyInstaller gsm_toolbox.spec --noconfirm

# step 2, about 4 minutes: wrap it into an installer
& "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe" installer\gsm_toolbox.iss
```

Results:

| Step | Produces |
|---|---|
| 1 | `dist\GSM_ToolBox\GSM_ToolBox.exe`, runnable straight away |
| 2 | `installer\Output\GSM_ToolBox_Setup_<version>.exe` |

To make the portable zip instead of an installer, compress the folder from step 1:

```powershell
Compress-Archive -Path dist\GSM_ToolBox -DestinationPath GSM_ToolBox-portable.zip
```

If you open the `.iss` in the Inno Setup Compiler window rather than running `ISCC.exe`,
use Build, then Compile. The same two-step order applies.

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
