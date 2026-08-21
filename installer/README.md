# Building the Windows installer

`gsm_toolbox.iss` is an [Inno Setup 6](https://jrsoftware.org/isinfo.php) script. It
packages a PyInstaller build that already exists. It does not create one.

Run the two steps in this order from the repository root:

```powershell
python -m PyInstaller gsm_toolbox.spec --noconfirm
& "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe" installer\gsm_toolbox.iss
```

The first step writes `dist\GSM_ToolBox`. The second reads that folder and writes
`installer\Output\GSM_ToolBox_Setup_<version>.exe`.

Compiling the script before the freeze step stops with a message telling you so. That is
the expected behaviour on a fresh clone, because a build directory of roughly 900 MB
cannot be kept in a git repository.

## What the installer does

It installs per user into `%LOCALAPPDATA%\Programs\GSM ToolBox`, so no administrator
rights are needed. It adds a Start menu entry, an optional desktop shortcut, and an
uninstaller. Downloaded databases in `%USERPROFILE%\.gsm_toolbox` are left alone when
uninstalling, because rebuilding them can take hours.

## Versioning

`AppVersion` near the top of the script must match `__version__` in
`gsm_toolbox/__init__.py`. The output filename carries the version, so an old installer
is never overwritten by a new build.
