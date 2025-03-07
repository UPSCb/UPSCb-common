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

# cmds container
cmds=()

# end of the boilerplate, logic goes below - instead of running cmds, add them to the cmds list
# e.g. 
cmds+=("echo Hello World
")

# end of logic, start of evalution. dry-run unless -d is provided on the cmdline

# shellcheck disable=SC2086
if [ ${DO} -eq 1 ]; then
    for j in $(seq 0 $((${#cmds[@]} - 1))); do
        eval "${cmds[$j]}"
    done
else
    echo "${cmds[@]}"
fi
