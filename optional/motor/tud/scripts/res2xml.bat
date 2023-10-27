::  RES2XML: Convert CFX5 result file to G+Smo parameterization
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
  if /i "%1" == "-e" set NUMELEVATE=%2 & shift
  if /i "%1" == "-o" set OUTPUTFILE=%2 & shift
  if /i "%1" == "-p" set PARAMFILE=%2 & shift
  if /i "%1" == "-r" set NUMREFINE=%2 & shift
  if /i "%1" == "-V" set VERBOSE=1
  shift
  :: Set input filename
  if not "%1" == "" (
    set INPUTFILE=%1
  goto getopts
  )

:: Set output filename
if .!OUTPUTFILE!==. set OUTPUTFILE=!INPUTFILE!

:: Set number of degree elevation steps
if .!NUMELEVATE!==. set NUMELEVATE=0

:: Set number of basis refinement steps
if .!NUMREFINE!==. set NUMREFINE=0

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
   if .!PARAMFILE!==. (
      echo !CMD! %CONVERTER% ^
          -i                  %INPUTFILE% ^
          -o                  %OUTPUTFILE% ^
          -e                  %NUMELEVATE ^
          -r                  %NUMREFINE ^
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
   ) else (
      echo !CMD! %CONVERTER% ^
          -i                  %INPUTFILE% ^
          -o                  %OUTPUTFILE% ^
          -p                  %PARAMFILE% ^
          -e                  %NUMELEVATE ^
          -r                  %NUMREFINE ^
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

if .!PARAMFILE!==. (
   !CMD! %CONVERTER% ^
    -i                  %INPUTFILE% ^
    -o                  %OUTPUTFILE% ^
    -e                  %NUMELEVATE ^
    -r                  %NUMREFINE ^
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
) else (
   !CMD! %CONVERTER% ^
    -i                  %INPUTFILE% ^
    -o                  %OUTPUTFILE% ^
    -p                  %PARAMFILE% ^
    -e                  %NUMELEVATE ^
    -r                  %NUMREFINE ^
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
echo Usage: res2xml.bat [OPTION]... FILE
echo Convert CFX5 result FILE to G+Smo XML parameterization.
echo Options:
echo   -h,-?             show this help
echo   -c configfile     use config file 'configfile' (default 'config.cfg')
echo   -d                run docker image (default no)
echo   -e numelevate     perform 'numelevate' degree elevation steps (default 0)
echo   -o outputfile     use output file 'outputfile'. If not given then the input file is overwritten
echo   -p paramfile      use parameter file 'paramfile'. If not given then the parameterization from the input file is adopted
echo   -r numrefine      perform 'numrefine' basis refinement steps (default 0)
echo   -V                verbose output (default no)
goto:eof
