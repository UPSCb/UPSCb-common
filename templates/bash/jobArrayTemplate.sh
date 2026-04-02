#!/bin/bash -l

# be safe
set -euo pipefail
#set -x

# variables
ARGs=0
DO=0
START=0
END=

# EDIT to adjust accordingly, also EDIT line 72
DATA=$(realpath ../../data)
PATTERN="*.txt"

# usage
export USAGETXT=\
"
Usage: $0 [options]

Purpose: The script is meant to ...

Options:
    -d do not just print, do
    -e end the array at a given task
    -h print this message
    -s start the array at task [default ${START}]

Note:
    The script is meant to be called as a job array.
"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts de:hs: option
do
    case "$option" in
        d) DO=1;;
		e) END=${OPTARG};;
        h) usage;;
		s) START=${OPTARG};;
		\?) ## unknown flag
		usage;;
        esac
done
shift $((OPTIND - 1))

# setup

# sanity
[[ ${ARGs} -eq 0 ]] && [[ $# -gt 0 ]] && abort "You provided arguments to a script that does not expect any"
[[ ${ARGs} -gt 0 ]] && [[ $# -ne ${ARGs} ]] && abort "This script expects ${ARGs} arguments"
[[ ${START} -lt 0 ]] && abort "-s cannot be smaller than 0"

# cmds container
cmds=()

# set the max number of tasks
read -r -a FILELIST <<< $(find "${DATA}" -mindepth 1 -maxdepth 1 -type f -name "${PATTERN}" -print0 | xargs -0)
[[ -z ${END} ]] && LEND=$((${#FILELIST[@]} - 1)) || LEND=${END}

# more sanity
[[ ${LEND} -lt ${START} ]] && abort "-e cannot be smaller than -s"
[[ ${LEND} -gt ${#FILELIST[@]} ]] && abort "-e cannot be larger than the maximum number of files to proceed: ${#FILELIST[@]}"

# job array submission
# EDIT THE FOLLOWING TO PASS THE RIGHT NUMBER OF ARGUMENTS AND ADJUST OPTIONS TO THE TASK SCRIPT
cmds+=("sbatch --array=${START}-${LEND} jobArrayTask.sh -d [options] ${DATA} <arguments>
")

# shellcheck disable=SC2086
if [ ${DO} -eq 1 ]; then
    for j in $(seq 0 $((${#cmds[@]} - 1))); do
        eval "${cmds[$j]}"
    done
else
    echo "${cmds[@]}"
	echo "
Dry-run is the default. Re-run with the -d option to submit the job.
"
fi
