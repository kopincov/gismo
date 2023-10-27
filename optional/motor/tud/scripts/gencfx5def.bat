::  GENCFX5DEF: Generate CFX5-solver definition
::
::  This file is part of the G+Smo library.
::
::  This Source Code Form is subject to the terms of the Mozilla Public
::  License, v. 2.0. If a copy of the MPL was not distributed with this
::  file, You can obtain one at http://mozilla.org/MPL/2.0/.
::
::  Author: M. Moller

@setlocal enableextensions enabledelayedexpansion
@echo off
set interactive=0

:: Set default configuration file
set CFG=%CD%\config.cfg

:: Get command line arguments
:getopts 
  if /i "%1" == "-h" goto Help
  if /i "%1" == "-?" goto Help
  if /i "%1" == "-c" set CFG=%2 & shift
  if /i "%1" == "-m" set MESHFILE=%2 & shift
  if /i "%1" == "-o" set OUTPUTFILE=%2 & shift
  if /i "%1" == "-V" set VERBOSE=1
  shift
  :: Set input filename
  if not "%1" == "" (
    set INPUTFILE=%1
  goto getopts
  )

:: Set output filename
if .!OUTPUTFILE!==. set OUTPUTFILE=!INPUTFILE!

:: Set mesh file
if .!MESHFILE!==. set MESHFILE=!~nINPUTFILE!.cfx5

:: Check that the config file exists
if not exist %CFG% (
  echo Config file %CFG% does not exist. Aborting.
  exit /b 1
)

:: Read config file
call:cfg_parser %CFG%

:: Import local variables from parser
:read_cfg
for /f "tokens=1* delims=@" %%a in ("%KEY_VAL_PAIR%") do (
  set %%a
  set KEY_VAL_PAIR=%%b
)
if not .!KEY_VAL_PAIR!==. goto read_cfg

:: Determine CFXROOT directory
if .!CFXROOT!==. (
  if not .!ANSYS190_DIR!==. (
    set CFXROOT=!ANSYS190_DIR!\..\CFX
  ) else if not .!ANSYS182_DIR!==. (
    set CFXROOT=!ANSYS182_DIR!\..\CFX
  ) else if not .!ANSYS181_DIR!==. (
    set CFXROOT=!ANSYS181_DIR!\..\CFX
  ) else if not .!ANSYS172_DIR!==. (
    set CFXROOT=!ANSYS172_DIR!\..\CFX
  ) else if not .!ANSYS171_DIR!==. (
    set CFXROOT=!ANSYS171_DIR!\..\CFX
  ) else if not .!ANSYS162_DIR!==. (
    set CFXROOT=!ANSYS162_DIR!\..\CFX
  ) else if not .!ANSYS161_DIR!==. (
    set CFXROOT=!ANSYS161_DIR!\..\CFX
  ) else (
    echo Unable to determine CFXROOT directory. Aborting.
    exit /b 1
  )
)

:: Check that executable cfx5info.exe is in PATH variable
if not exist !CFXROOT!\bin\cfx5info.exe (
  echo Command cfx5info.exe is not in PATH variable. Aborting.
  exit /b 1
)

:: Determine CFX5 version
for /f %%a in ('"!CFXROOT!\bin\cfx5info.exe" -version') do set CFXVERSION=%%a

:: Determine absolute paths
for /f %%a in ("!INPUTFILE!") do set INPUTFILE=%%~fa
for /f %%a in ("!MESHFILE!") do set MESHFILE=%%~fa
for /f %%a in ("!OUTPUTFILE!") do set OUTPUTFILE=%%~fa

:: Generate CFX5 session file
set TMPFILE=tmp-%RANDOM%.pre
@echo COMMAND FILE: > !TMPFILE!
@echo   CFX Pre Version = !CFXVERSION! >> !TMPFILE!
@echo END >> !TMPFILE!
@echo ^>load filename=!INPUTFILE!, ^\ >> !TMPFILE!
@echo mode=def, recoverSession=no, replaceFlow=yes, overwrite=yes >> !TMPFILE!
@echo ^>gtmAction, op=deleteAssembly, Assembly >> !TMPFILE!
@echo ^>gtmImport filename=!MESHFILE!, ^\ >> !TMPFILE!
@echo type=Generic, units=!cfx5pre_UNIT!, genOpt= ^-n, ^\ >> !TMPFILE!
@echo nameStrategy= Assembly >> !TMPFILE!
@echo ^>writeCaseFile filename=!OUTPUTFILE!, ^\ >> !TMPFILE!
@echo operation=start solver interactive >> !TMPFILE!

:: Check that executable cfx5pre.exe is in PATH variable
if not exist !CFXROOT!\bin\cfx5pre.exe (
  echo Command cfx5pre.exe is not in PATH variable. Aborting.
  exit /b 1
)

:: Call CFX5-pre in batch mode
if /i "!VERBOSE!"=="1" (
  echo "Running..."
  echo "!CFXROOT!\bin\cfx5pre" -batch !TMPFILE!
)
"!CFXROOT!\bin\cfx5pre" -batch !TMPFILE!

:: Delete CFX5 session file
del !TMPFILE!

if "%interactive%"=="0" pause
endlocal
exit /b 0
goto:eof

:: Parser for configuration file (Windows INI format)
:cfg_parser
@setlocal enableextensions enabledelayedexpansion
@echo off
set file=%~1
set section=[%~2]
set key=%~3
set currsection=

set end=
for /f "usebackq delims=" %%a in ("!file!") do (
    set ln=%%a
    if "x!ln:~0,1!"=="x[" (
        set currsection=!ln!
    ) else (
        for /f "tokens=1,2 delims==" %%b in ("!ln!") do (
            set currkey=%%b
            set currval=%%c
      
      for /f "tokens=* delims= " %%a in ("!currkey!") do set currkey=%%a
      for /f "tokens=1 delims= " %%a in ("!currkey!") do set currkey=%%a
      for /f "tokens=* delims= " %%a in ("!currval!") do set currval=%%a
      
      set KEY_VAL_PAIR=!currsection:~1,-1!_!currkey!=!currval!@!KEY_VAL_PAIR!
        )
    )
)

for /f "tokens=*" %%a in ("!KEY_VAL_PAIR!") do (
  endlocal & set "KEY_VAL_PAIR=%%~a"
)
goto:eof

:: Print help screen
: Help
echo Usage: cfx5updatemesh.bat [OPTION]... FILE
echo Generate CFX5-solver definition FILE.
echo Options:
echo   -h,-?             show this help
echo   -c configfile     use config file 'configfile' (default 'config.cfg')
echo   -m meshfile       use mesh file 'meshfile' (default 'mesh.cfx5')
echo   -o outputfile     use output file 'outputfile'. If not given then the input file is overwritten
echo   -V                verbose output (default no)
goto:eof
