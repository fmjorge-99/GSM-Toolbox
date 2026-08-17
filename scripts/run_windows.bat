@echo off
REM Launch GSM ToolBox from source on Windows, creating the environment on first run.
setlocal
set "HERE=%~dp0.."
set "VENV=%HERE%\.venv"

if not exist "%VENV%\Scripts\python.exe" (
    echo Creating virtual environment...
    py -3.12 -m venv "%VENV%" 2>nul || py -3.11 -m venv "%VENV%" 2>nul || python -m venv "%VENV%"
    if errorlevel 1 (
        echo error: could not create a virtual environment. Python 3.10-3.12 required.
        exit /b 1
    )
    "%VENV%\Scripts\python.exe" -m pip install --upgrade pip
    echo Installing dependencies ^(a few minutes the first time^)...
    "%VENV%\Scripts\python.exe" -m pip install -r "%HERE%\requirements.txt"
)

"%VENV%\Scripts\python.exe" "%HERE%\run_gsm_toolbox.py" %*
endlocal
