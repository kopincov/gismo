#!/bin/bash

#  GENCFX5DEF: Generate CFX5-solver definition
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
while getopts "h?c:m:o:V" opt; do
    case "$opt" in
        h|\?)
            echo "Usage: gencfx5def.sh [OPTION]... FILE";
            echo "Generate CFX5-solver definition FILE."
            echo "Options:"
            echo "  -h,-?             show this help"
            echo "  -c configfile     use config file 'configfile' (default 'config.cfg')"
            echo "  -m meshfile       use mesh file 'meshfile' (default 'mesh.cfx5')"
            echo "  -o outputfile     use output file 'outputfile'. If not given then the input file is overwritten."
            echo "  -V                verbose output (default no)"
            exit 0
            ;;
        c)  CFG=$OPTARG
            ;;
        m)  MESHFILE=$OPTARG
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
    OUTPUTFILE=$INPUTFILE
fi

# Set mesh filename
if [ -z "${MESHFILE+x}" ]; then
    MESHFILE=$(echo $INPUTFILE | cut -f 1 -d '.').cfx5
fi

# Check that the config file exists
if [ ! -e "$CFG" ]; then
    echo "Config file '$CFG' does not exist. Aborting."
    exit 1;
fi

# Read config file
cfg_parser $CFG

# Determine CFXROOT directory
if [ -z "${CFXROOT+x}" ]; then
    if [ -n "$ANSYS190_DIR" ]; then
        CFXROOT=$ANSYS190_DIR/../CFX
    elif [ -n "$ANSYS182_DIR" ]; then
        CFXROOT=$ANSYS182_DIR/../CFX
    elif [ -n "$ANSYS181_DIR" ]; then
        CFXROOT=$ANSYS181_DIR/../CFX
    elif [ -n "$ANSYS172_DIR" ]; then
        CFXROOT=$ANSYS172_DIR/../CFX
    elif [ -n "$ANSYS171_DIR" ]; then
        CFXROOT=$ANSYS171_DIR/../CFX
    elif [ -n "$ANSYS162_DIR" ]; then
        CFXROOT=$ANSYS162_DIR/../CFX
    elif [ -n "$ANSYS161_DIR" ]; then
        CFXROOT=$ANSYS161_DIR/../CFX
    else
        echo "Unable to determine CFXROOT directory. Aborting."
        exit 1;
    fi
fi

# Check that cfx5info executable is in PATH variable
command -v $CFXROOT/bin/cfx5info >/dev/null 2>&1 || \
    { echo "Command cfx5info is not in PATH variable. Aborting." >&2; exit 1; }

# Determine CFX5 version
CFXVERSION=`$CFXROOT/bin/cfx5info -version`

# Determine absolute paths
INPUTFILE=$(realpath "$INPUTFILE")
MESHFILE=$(realpath "$MESHFILE")
OUTPUTFILE=$(realpath "$OUTPUTFILE")

# Generate CFD5 session file
TMPFILE=`mktemp tmp-XXXXXXXX.pre`
echo "COMMAND FILE:" >> $TMPFILE
echo "  CFX Pre Version = $CFXVERSION" >> $TMPFILE
echo "END" >> $TMPFILE
echo ">load filename=$INPUTFILE, \\" >> $TMPFILE
echo "mode=def, recoverSession=no, replaceFlow=yes, overwrite=yes" >> $TMPFILE
echo ">gtmAction, op=deleteAssembly, Assembly" >> $TMPFILE
echo ">gtmImport filename=$MESHFILE, \\" >> $TMPFILE
echo "type=Generic, units=$cfx5pre_UNIT, genOpt= -n, \\" >> $TMPFILE
echo "nameStrategy= Assembly" >> $TMPFILE
echo ">writeCaseFile filename=$OUTPUTFILE, \\" >> $TMPFILE
echo "operation=start solver interactive" >> $TMPFILE

# Check that cfx5pre executable is in PATH variable
command -v $CFXROOT/bin/cfx5pre >/dev/null 2>&1 || \
    { echo "Command cfx5pre is not in PATH variable. Aborting." >&2; exit 1; }

# Call CFX5-pre in batch mode
if [ -z "${VERBOSE+x}" ]; then
    command $CFXROOT/bin/cfx5pre -batch $TMPFILE
else
    echo "Running..."
    echo $CFXROOT/bin/cfx5pre -batch $TMPFILE
    command $CFXROOT/bin/cfx5pre -batch $TMPFILE
fi

# Delete CFX5 session file
rm -f $TMPFILE

exit 0;
