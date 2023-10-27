#!/bin/bash

#  XML2CFX5: Convert G+Smo parameterization to CFX5 mesh
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

# A POSIX variable
OPTIND=1
while getopts "h?c:do:V" opt; do
    case "$opt" in
        h|\?)
            echo "Usage: xml2cfx5.sh [OPTION]... FILE";
            echo "Convert G+Smo XML parameterization FILE into CFX5 mesh."
            echo "Options:"
            echo "  -h,-?             show this help"
            echo "  -c configfile     use config file 'configfile' (default 'config.cfg')"
            echo "  -d                run docker image (default no)"
            echo "  -o outputfile     use output file 'outputfile'. If not given then the CFX5 mesh API is used"
            echo "  -V                verbose output (default no)"
            exit 0
            ;;
        c)  CFG=$OPTARG
            ;;
        d)  DOCKER=1
            ;;
        o)  OUTPUTFILE=$OPTARG
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
    OUTPUTFILE=-
    CONVERTER=--xml_cfx5
else
    CONVERTER=--xml_icem
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
    command $CMD $CONVERTER \
        -i                  $INPUTFILE \
        -o                  $OUTPUTFILE \
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
    echo "Running..."
    echo $CMD $CONVERTER \
        -i                  $INPUTFILE \
        -o                  $OUTPUTFILE \
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
    command $CMD $CONVERTER \
        -i                  $INPUTFILE \
        -o                  $OUTPUTFILE \
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
exit 0;
