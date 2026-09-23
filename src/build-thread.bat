@echo off
setlocal

set SRC=src\solution-thread.cpp
set OUT=solution-thread.exe

g++ "%SRC%" -O2 -o "%OUT%" -lgmp
if errorlevel 1 (
    echo Build failed.
    exit /b 1
)

echo Built %OUT%

if not "%~1"=="" (
    "%OUT%" %*
)

endlocal
