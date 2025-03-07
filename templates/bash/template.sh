#!/bin/bash

# failsafe
set -eu -o pipefail
#set -x

# variables
ARGs=0
DO=0

# usage
export USAGETXT=\
"
Usage: $0 [options] <arguments>

Purpose: The script ...

Options:
    -d do not just print, do
    -h print this message
    
Note:
    The script ...
"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts dh option
do
        case "$option" in
        d) DO=1;;
        h) usage;;
		\?) ## unknown flag
		usage;;
        esac
done
shift $((OPTIND - 1))

# setup

# sanity
[[ ${ARGs} -gt 0 ]] && [[ $# -ne ${ARGs} ]] && abort "This script expects ${ARGs} arguments"

# end of the boilerplate, logic goes below
