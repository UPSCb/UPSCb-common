#!/bin/bash -l
# EDIT to add default SBATCH (and remove this line)
#SBATCH -A ...
#SBATCH -c ...
#SBATCH -t ...

# be safe
set -euo pipefail
#set -x

# variables
ARGs=2
DO=0
SEQUENTIAL=0

# EDIT to adjust accordingly, also EDIT line 72
PATTERN="*.txt"

# usage
export USAGETXT=\
"
Usage: $0 [options] <input directory> <arguments>

Purpose: The script is meant to ...

Options:
    -d do not just print, do
    -h print this message
    -s sequential use, input becomes files instead of directories

Note: It is meant to run as a job array, use the -s option above to run it on a single file.
"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts dhs option
do
    case "$option" in
        d) DO=1;;
        h) usage;;
		s) SEQUENTIAL=1;;
		\?) ## unknown flag
			usage;;
    esac
done
shift $((OPTIND - 1))

# sanity
[[ ${ARGs} -gt 0 ]] && [[ $# -ne ${ARGs} ]] && abort "This script expects ${ARGs} arguments"

[[ ${SEQUENTIAL} -eq 0 ]] && [[ ! -d "${1}" ]] && abort "The first argument needs to be an existing directory"

[[ ${SEQUENTIAL} -eq 1 ]] && [[ ! -f "${1}" ]] && abort "The first argument needs to be an existing file"

# setup
IN="${1}" && shift

# NOTE - if your files are named <prefix>.<task>, you can comment this block and use the single line commented out after the block
if [ ${SEQUENTIAL} -eq 0 ]; then
	read -r -a LIST <<< $(find "${IN}" -mindepth 1 -maxdepth 1 -type f -name "${PATTERN}" -print0 | xargs -0 | sort)
	# the :-0 syntax is just to allow for testing without SLURM, setting the default index value to 0 (i.e. job array task 0)
	IN=${LIST[${SLURM_ARRAY_TASK_ID:-0}]}
fi
# HERE - the alternative to the block above relies on PATTERN to be set to the prefix
# [[ ${SEQUENTIAL} -eq 0 ]] && IN=${IN}/${PATTERN}.${SLURM_ARRAY_TASK_ID:-0}


# cmds container
cmds=()

# end of the boilerplate, logic goes below - instead of running cmds, add them to the cmds list
# e.g. 
cmds+=("echo Hello World, I'm ${IN}
")

# shellcheck disable=SC2086
if [ ${DO} -eq 1 ]; then
    for j in $(seq 0 $((${#cmds[@]} - 1))); do
        eval "${cmds[$j]}"
    done
else
    echo "${cmds[@]}"
fi
