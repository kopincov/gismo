::  XML2CFX5: Convert G+Smo parameterization to CFX5 mesh
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
  if /i "%1" == "-d" set DOCKER=1
  if /i "%1" == "-o" set OUTPUTFILE=%2 & shift
  if /i "%1" == "-V" set VERBOSE=1
  shift
  :: Set input filename
  if not "%1" == "" (
    set INPUTFILE=%1
  goto getopts
  )

:: Set output filename
if .!OUTPUTFILE!==. (
   set OUTPUTFILE=-
   set CONVERTER=--xml_cfx5
) else (
  set CONVERTER=--xml_icem
)

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

:: Replace $(pwd) by %cd% (if required)
set docker_CMD=!docker_CMD:$(pwd)=%cd%!

:: Check that executable is in PATH variable
if /i "!DOCKER!"=="1" (
    where /q docker 
  if ERRORLEVEL 1 (
      echo Docker executable is not in PATH variable. Aborting.
      echo /b 1
    )
    set CMD=!docker_CMD! !docker_IMAGE! !gismo_CMD!
    ) else (
    where /q "!gismo_CMD!" 
  if ERRORLEVEL 1 (
      echo Command !gismo_CMD! is not in PATH variable. Aborting.
      echo /b 1
    )
    set CMD=!gismo_CMD!
)

:: Run executable with parameters from config file and command line
if /i "!VERBOSE!"=="1" (
  echo "Running..."
  echo !CMD! %CONVERTER% ^
    -i                  %INPUTFILE% ^
    -o                  %OUTPUTFILE% ^
    --nleft             !gismo_NLEFT! ^
    --nmiddle           !gismo_NMIDDLE! ^
    --nright            !gismo_NRIGHT! ^
    --ntop              !gismo_NTOP! ^
    --nbottom           !gismo_NBOTTOM! ^
    --nz                !gismo_NZ! ^
    --xscale            !gismo_XSCALE! ^
    --yscale            !gismo_YSCALE! ^
    --zscale            !gismo_ZSCALE! ^
    --xleftbottom       "!gismo_XLEFTBOTTOM!" ^
    --xmiddlebottom     "!gismo_XMIDDLEBOTTOM!" ^
    --xrightbottom      "!gismo_XRIGHTBOTTOM!" ^
    --xleftmiddle       "!gismo_XLEFTMIDDLE!" ^
    --xrightmiddle      "!gismo_XRIGHTMIDDLE!" ^
    --xlefttop          "!gismo_XLEFTTOP!" ^
    --xmiddletop        "!gismo_XMIDDLETOP!" ^
    --xrighttop         "!gismo_XRIGHTTOP!" ^
    --yleftbottomleft   "!gismo_YLEFTBOTTOMLEFT!" ^
    --yleftbottomright  "!gismo_YLEFTBOTTOMRIGHT!" ^
    --ylefttopleft      "!gismo_YLEFTTOPLEFT!" ^
    --ylefttopright     "!gismo_YLEFTTOPRIGHT!" ^
    --yrighttopleft     "!gismo_YRIGHTTOPLEFT!" ^
    --yrighttopright    "!gismo_YRIGHTTOPRIGHT!" ^
    --yrightbottomright "!gismo_YRIGHTBOTTOMRIGHT!" ^
    --yrightbottomleft  "!gismo_YRIGHTBOTTOMLEFT!"
)
!CMD! %CONVERTER% ^
    -i                  %INPUTFILE% ^
    -o                  %OUTPUTFILE% ^
    --nleft             !gismo_NLEFT! ^
    --nmiddle           !gismo_NMIDDLE! ^
    --nright            !gismo_NRIGHT! ^
    --ntop              !gismo_NTOP! ^
    --nbottom           !gismo_NBOTTOM! ^
    --nz                !gismo_NZ! ^
    --xscale            !gismo_XSCALE! ^
    --yscale            !gismo_YSCALE! ^
    --zscale            !gismo_ZSCALE! ^
    --xleftbottom       "!gismo_XLEFTBOTTOM!" ^
    --xmiddlebottom     "!gismo_XMIDDLEBOTTOM!" ^
    --xrightbottom      "!gismo_XRIGHTBOTTOM!" ^
    --xleftmiddle       "!gismo_XLEFTMIDDLE!" ^
    --xrightmiddle      "!gismo_XRIGHTMIDDLE!" ^
    --xlefttop          "!gismo_XLEFTTOP!" ^
    --xmiddletop        "!gismo_XMIDDLETOP!" ^
    --xrighttop         "!gismo_XRIGHTTOP!" ^
    --yleftbottomleft   "!gismo_YLEFTBOTTOMLEFT!" ^
    --yleftbottomright  "!gismo_YLEFTBOTTOMRIGHT!" ^
    --ylefttopleft      "!gismo_YLEFTTOPLEFT!" ^
    --ylefttopright     "!gismo_YLEFTTOPRIGHT!" ^
    --yrighttopleft     "!gismo_YRIGHTTOPLEFT!" ^
    --yrighttopright    "!gismo_YRIGHTTOPRIGHT!" ^
    --yrightbottomright "!gismo_YRIGHTBOTTOMRIGHT!" ^
    --yrightbottomleft  "!gismo_YRIGHTBOTTOMLEFT!"

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
echo Usage: xml2cfx5.bat [OPTION]... FILE
echo Convert G+Smo XML parameterization FILE into CFX5 mesh.
echo Options:
echo   -h,-?             show this help
echo   -c configfile     use config file 'configfile' (default 'config.cfg')
echo   -d                run docker image (default no)
echo   -o outputfile     use output file 'outputfile'. If not given then the CFX5 mesh API is used
echo   -V                verbose output (default no)
goto:eof
