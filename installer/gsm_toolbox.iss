; Inno Setup script for GSM ToolBox.
; Produces a single compact Setup .exe that installs the one-folder PyInstaller
; build (dist\GSM_ToolBox) plus all its bundled dependencies (Python runtime,
; COBRApy, RDKit, SCIP, the offline universal database, …) on Windows 10+.
;
; This script packages an existing PyInstaller build. It does NOT create one, so the
; application must be frozen first or there is nothing to install:
;
;   1.  pip install -r requirements.txt pyinstaller
;   2.  python -m PyInstaller gsm_toolbox.spec --noconfirm
;   3.  "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" installer\gsm_toolbox.iss
;
; Step 2 writes dist\GSM_ToolBox, which step 3 reads. Output lands in
; installer\Output\GSM_ToolBox_Setup_<version>.exe.
;
; The app is self-contained, so no admin rights are required: it installs per-user
; into %LOCALAPPDATA%\Programs by default. LZMA2/ultra compression keeps the
; installer as small as possible.

#define AppName "GSM ToolBox"
#define AppVersion "0.3.10"
#define AppPublisher "Jorge Fernandez Mendez"
#define AppExe "GSM_ToolBox.exe"

[Setup]
AppId={{7F3C1A2E-6B4D-4E9A-9C21-A1B2C3D4E5F6}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
DefaultDirName={autopf}\GSM ToolBox
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
; Per-user install so no admin prompt is needed; the bundle is fully self-contained.
PrivilegesRequiredOverridesAllowed=dialog commandline
PrivilegesRequired=lowest
OutputDir=Output
OutputBaseFilename=GSM_ToolBox_Setup_{#AppVersion}
SetupIconFile=..\gsm_toolbox\resources\icons\app_icon.ico
UninstallDisplayIcon={app}\{#AppExe}
WizardStyle=modern
Compression=lzma2/ultra64
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional icons:"

; Fail with an explanation rather than "No files found matching". A fresh clone has no
; dist folder, because a 900 MB build directory cannot live in a git repository, and the
; bare Inno error sends people looking for a missing file instead of a missing step.
#define FrozenDir "..\dist\GSM_ToolBox"
#if !DirExists(AddBackslash(SourcePath) + FrozenDir)
  #error The frozen application was not found at dist\GSM_ToolBox. Run "python -m PyInstaller gsm_toolbox.spec --noconfirm" from the repository root first, then compile this script again. See INSTALL.md, section 4.
#endif

[Files]
; Recursively bundle the entire one-folder PyInstaller build.
Source: "{#FrozenDir}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "Launch {#AppName}"; Flags: nowait postinstall skipifsilent
