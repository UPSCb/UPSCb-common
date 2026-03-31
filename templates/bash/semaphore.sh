#!/bin/bash -l

# failsafe
set -eu -o pipefail
#set -x

# variables
CLEAN=0
DO=0
FORCE=0
RESUME=0
WAIT=10

# EDIT these
POOLSIZE=2
WORKDIR=""

# usage
export USAGETXT=\
"
Usage: $0 [options]

Purpose: The script will ...

Options:
    -c clean the pool, especially after a dry-run
    -d do not just print, do
    -f force rerunning
    -h print this message
    -r resume from a previous run (i.e. the pool has already been used and populated)
    -w wait time for the lock (default: ${WAIT} seconds)

Note:
    The script uses a semaphore to manage ...
    
    -c and -r are mutually exclusive
"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts cdfhl:p:rw: option
do
    case "$option" in
        c) CLEAN=1;;
        d) DO=1;;
        f) FORCE=1;;
        h) usage;;
        r) RESUME=1;;
        w) WAIT=${OPTARG};;
		\?) ## unknown flag
		usage;;
    esac
done
shift $((OPTIND - 1))

# sanity
# EDIT to add checks

# derivatives
LOG="${WORKDIR}/log"
POOL="${WORKDIR}/pool"
TMP="${WORKDIR}/tmp"

# setup
[[ ! -d "${LOG}" ]] && mkdir -p "${LOG}"
[[ ! -d "${POOL}" ]] && mkdir -p "${POOL}"
[[ ! -d "${TMP}" ]] && mkdir -p "${TMP}"

POOLFILE="${POOL}/size.txt"

if [ ! -f "$POOLFILE" ]; then
	echo "${POOLSIZE}" > "${POOLFILE}"
else
	# sanity
	[[ ${CLEAN} -eq 1 ]] && [[ ${RESUME} -eq 1 ]] && abort "-c and -r are mutually exclusive"
	# clean
	[[ ${CLEAN} -eq 1 ]] && echo "${POOLSIZE}" > "${POOLFILE}"
	# resume
	[[ ${CLEAN} -eq 0 ]] && [[ ${RESUME} -eq 0 ]] && abort "The pool is already in use. Either clean (-c) or resume (-r)"
fi        

# preparation
# WRITE and EDIT so you get a list object having the list of files to process (or project, whatever...)
for l in $(seq 0 $((${#list[@]} - 1))); do

	ls=${list[$l]}

	# skip if done
	if [ -f "${LOG}/$(basename "${ls}").completed" ] && [ ${FORCE} -eq 0 ]; then
		echo "Skipping ${ls} as it is already processed. Use \"-l ${ls} -f\" to force redo it."
	else

		while true; do

			tmp=$(mktemp -p "${TMP}")

			# create a lock (flock) and read the size of the POOL
			# from https://unix.stackexchange.com/questions/70/what-unix-commands-can-be-used-as-a-semaphore-lock
			# and man flock
			(
				# Wait for lock (assigned to file descriptor 200)
				# See https://stackoverflow.com/questions/5256599/what-are-file-descriptors-explained-in-simple-terms for explanation of file descriptors
				flock -x -w "${WAIT}" 200

				# read the file content
				read -r POOLSIZE < "${POOLFILE}"

				[[ ${POOLSIZE} -gt 0 ]] && echo 1 > "${tmp}" || echo 0 > "${tmp}"
			
			) 200> "${POOLFILE/.txt/.lock}"

			# read the result of the lock
			read -r RESULT < "${tmp}"
			rm "${tmp}"

			# shellcheck disable=SC2086
			if [ ${RESULT} -eq 0 ]; then
				# if it is below 0, wait
				sleep 60
				continue
			else
				# otherwise, submit the job
				(bash process.sh -f "${FORCE}" -d "${DO}" -w "${WAIT}" "${WORKDIR}" "${ls}")&
				(
					# update the pool size
					flock -x -w "${WAIT}" 200
					read -r POOLSIZE < "${POOLFILE}"
					echo $((POOLSIZE - 1)) > "${POOLFILE}"
				) 200> "${POOLFILE/.txt/.lock}"
				break
			fi
		done
	fi
done
