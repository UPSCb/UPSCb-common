#!/bin/bash -l

# failsafe
set -eu -o pipefail
#set -x

# variables
ARGs=2
DO=0
FORCE=0
WAIT=10

# usage
export USAGETXT=\
"
Usage: $0 [options] <working directory> <file list>

Purpose: The script will create ...
Options:
    -d do not just print, do
    -f force redo
    -w wait time for the lock (default: ${WAIT} seconds)
"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts d:f:w: option
do
        case "$option" in
        d) DO=${OPTARG};;
        f) FORCE=${OPTARG};;
        w) WAIT=${OPTARG};;
		\?) ## unknown flag
		usage;;
        esac
done
shift $((OPTIND - 1))

# sanity and arguments
[[ ${ARGs} -eq 0 ]] && [[ $# -gt 0]] && abort "You provided arguments to a script that does not expect any"
[[ ${ARGs} -gt 0 ]] && [[ $# -ne ${ARGs} ]] && abort "This script expects ${ARGs} arguments"


[[ ${DO} -ne 0 ]] && [[ ${DO} -ne 1 ]] && abort "-d should have a value of 0 or 1"

[[ ${FORCE} -ne 0 ]] && [[ ${FORCE} -ne 1 ]] && abort "-f should have a value of 0 or 1"

[[ ! -d $1 ]] && abort "The working directory $1 does not exist"
WORKDIR=$1
shift

[[ ! -f $1 ]] && abort "The file $1 does not exist"
FILE=$1
shift

# derivatives
LOG="${WORKDIR}/log"
POOL="${WORKDIR}/pool"
TMP="${WORKDIR}/tmp"
POOLFILE="${POOL}/size.txt"

# setup
[[ ! -d "${TMP}" ]] && mkdir -p "${TMP}"

# end of first boiler plate code section
# WRITE and EDIT to add yout jobs / commands to run
# I like to call separate task scripts
# first task, e.g. tar
[[ FORCE -eq 1 ]] || [[ ! -f "${ARCHIVE}" ]] && bash <COMMAND.sh> [options] <arguments>

# second task, e.g. md5sum
[[ FORCE -eq 1 ]] || [[ ! -f "${LOG}/${NAME}.completed" ]] && bash <COMMAND.sh> [options] <arguments>

# third task, e.g. rsync
[[ FORCE -eq 1 ]] || [[ ! -f "${LOG}/${NAME}.completed" ]] && bash <COMMAND.sh> [options] <arguments>

# mark as complete
# EDIT / create a NAME variable to set the completion flag
touch "$LOG/${NAME}.completed"


# second boiler plate section
# update the poolsize
# shellcheck disable=SC2094
(
    flock -x -w "${WAIT}" 200

    # read the file content
    read -r poolsize < "${POOLFILE}"

    # update it
    POOLSIZE=$((poolsize + 1))
    echo ${POOLSIZE} > "${POOLFILE}"

) 200> "${POOLFILE/.txt/.lock}"

# end of the second boiler plate section
# EDIT to add code to clean up / tear down
# e.g. rm "${TMP}"
