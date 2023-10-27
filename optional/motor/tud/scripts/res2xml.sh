#!/bin/bash

#  RES2XML: Convert CFX5 result file to G+Smo parameterization
#
#  This file is part of the G+Smo library.
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Author: M. Moller

# Parser for configuration file (Windows INI format)
function cfg_parser() {
    set +a
    while read p; do
        reSec='^\[(.*)\]$'
        #reNV='[ ]*([^ ]*)+[ ]*=(.*)'     #Remove only spaces around name
        reNV='[ ]*([^ ]*)+[ ]*=[ ]*(.*)'  #Remove spaces around name and spaces before value
        if [[ $p =~ $reSec ]]; then
            section=${BASH_REMATCH[1]}
        elif [[ $p =~ $reNV ]]; then
            sNm=${section}_${BASH_REMATCH[1]}
            sVa=${BASH_REMATCH[2]}
            set -a
            eval "$(echo "$sNm"=\""$sVa"\")"
            set +a
        fi
    done < $1
}

# Set default configuration file
CFG="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/config.cfg

# A POSIX variable<
OPTIND=1
while getopts "h?c:de:o:p:r:V" opt; do
    case "$opt" in
        h|\?)
            echo "Usage: res2xml.sh [OPTION]... FILE";
            echo "Convert CFX5 result FILE to G+Smo XML parameterization."
            echo "Options:"
            echo "  -h,-?             show this help"
            echo "  -c configfile     use config file 'configfile' (default 'config.cfg')"
            echo "  -d                run docker image (default no)"
            echo "  -e numelevate     perform 'numelevate' degree elevation steps (default 0)"
            echo "  -o outputfile     use output file 'outputfile'. If not given then the input file is overwritten"
            echo "  -p paramfile      use parameter file 'paramfile'. If not given then the parameterization from the input file is adopted"
            echo "  -r numrefine      perform 'numrefine' basis refinement steps (default 0)"
            echo "  -V                verbose output (default no)"
            exit 0
            ;;
        c)  CFG=$OPTARG
            ;;
        d)  DOCKER=1
            ;;
        e)  NUMELEVATE=$OPTARG
            ;;
        o)  OUTPUTFILE=$OPTARG
            ;;
        p)  PARAMFILE=$OPTARG
            ;;
        r)  NUMREFINE=$OPTARG
            ;;
        V)  VERBOSE=1
            ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# Set input filename
INPUTFILE=$1

# Set output filename
if [ -z "${OUTPUTFILE+x}" ]; then
    OUTPUTFILE=$INPUTFILE
fi

# Set number of degree elevation steps
if [ -z "${NUMELEVATE+x}" ]; then
    NUMELEVATE=0
fi

# Set number of basis refinement steps
if [ -z "${NUMREFINE+x}" ]; then
    NUMREFINE=0
fi

# Check that the config file exists
if [ ! -e "$CFG" ]; then
    echo "Config file '$CFG' does not exist. Aborting."
    exit 1;
fi

# Read config file
cfg_parser $CFG

# Check that executable is set
if [ -z "${gismo_CMD+x}" ]; then
    echo "Variable \$gismo_CMD is not set in config file '$CFG'. Aborting."
    exit 1;
fi

# Check that executable is in PATH variable
if [ -z "${DOCKER+x}" ]; then
    command -v $gismo_CMD >/dev/null 2>&1 || \
        { echo "Command $gismo_CMD is not in PATH variable. Aborting." >&2; exit 1; }
    CMD=$gismo_CMD;
else
    command -v $DOCKER_CMD >/dev/null 2>&1 || \
        { echo "Docker executable is not in PATH variable. Aborting." >&2; exit 1; }
    CMD="$docker_CMD $docker_IMAGE $gismo_CMD"
fi

# Run executable with parameters from config file and command line
if [ -z "${VERBOSE+x}" ]; then
    if [ -z "${PARAMFILE+x}" ]; then
        command $CMD --res_xml \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
        > /dev/null 2>&1 || { echo "An error occured. Aborting." >&2; exit 1; }
    else
        command $CMD --res_xml \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -p                  $PARAMFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
        > /dev/null 2>&1 || { echo "An error occured. Aborting." >&2; exit 1; }
    fi
else
    echo "Running..."
    if [ -z "${PARAMFILE+x}" ]; then
        echo  $CMD --res_xml  \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
        command $CMD --res_xml \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
    else
        echo  $CMD --res_xml  \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -p                  $PARAMFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
        command $CMD --res_xml \
            -i                  $INPUTFILE \
            -o                  $OUTPUTFILE \
            -p                  $PARAMFILE \
            -e                  $NUMELEVATE \
            -r                  $NUMREFINE \
            --nleft             "$gismo_NLEFT" \
            --nmiddle           "$gismo_NMIDDLE" \
            --nright            "$gismo_NRIGHT" \
            --ntop              "$gismo_NTOP" \
            --nbottom           "$gismo_NBOTTOM" \
            --nz                "$gismo_NZ" \
            --xscale            "$gismo_XSCALE" \
            --yscale            "$gismo_YSCALE" \
            --zscale            "$gismo_ZSCALE" \
            --xleftbottom       "$gismo_XLEFTBOTTOM" \
            --xmiddlebottom     "$gismo_XMIDDLEBOTTOM" \
            --xrightbottom      "$gismo_XRIGHTBOTTOM" \
            --xleftmiddle       "$gismo_XLEFTMIDDLE" \
            --xrightmiddle      "$gismo_XRIGHTMIDDLE" \
            --xlefttop          "$gismo_XLEFTTOP" \
            --xmiddletop        "$gismo_XMIDDLETOP" \
            --xrighttop         "$gismo_XRIGHTTOP" \
            --yleftbottomleft   "$gismo_YLEFTBOTTOMLEFT" \
            --yleftbottomright  "$gismo_YLEFTBOTTOMRIGHT" \
            --ylefttopleft      "$gismo_YLEFTTOPLEFT" \
            --ylefttopright     "$gismo_YLEFTTOPRIGHT" \
            --yrighttopleft     "$gismo_YRIGHTTOPLEFT" \
            --yrighttopright    "$gismo_YRIGHTTOPRIGHT" \
            --yrightbottomright "$gismo_YRIGHTBOTTOMRIGHT" \
            --yrightbottomleft  "$gismo_YRIGHTBOTTOMLEFT"
    fi
fi
exit 0;
