#!/bin/bash
#SBATCH -p main
#SBATCH -c 1
#SBATCH --mem=16GB
#SBATCH -t 2-00:00:00

# failsafe
set -eu -o pipefail
#set -x

# variables
ARGs=1
BASE=0
CONTAINER="$(realpath "../singularity/R-4.4.3.sif")"
DO=0
RLIB="/mnt/picea/Modules/apps/compilers/R/4.4.3/lib/R/library"


# usage
export USAGETXT=\
"
Usage: $0 [options] <arguments>

Purpose: The script will use Rscript to run the provided R script

Options:
    -b use the base R installation only. Mutually exclusive with -l. Default is off.
    -c use a different R container, defaults to ${CONTAINER}. Provide the path to the container.
    -d do not just print, do
    -h print this message
    -l provide the R package library. For using only base R, set -b. Default is ${RLIB}
    
Note:
    * The script expects to have a directory in the parent directory of this script called _singularity_ that contains the singularity container to use. Ideally you would copy this script from the UPSCb-common template collection to your project pipeline directory.
    * The script similarly expects the UPSCb-common helper file found in UPSCb-common/src/bash/functions.sh to be copied in the same directory as this script.
    * If your R script expects arguments pass them as extra arguments to this script and make sure to update the ARGS variable accordingly line 8"

# helper function
# shellcheck disable=SC1091
source functions.sh

# handle the options
while getopts bc:dhl: option
do
        case "$option" in
        b) BASE=1
          RLIB="";;
        c) CONTAINER=${OPTARG};;
        d) DO=1;;
        h) usage;;
        l) RLIB=${OPTARG};;
		    \?) ## unknown flag
		usage;;
        esac
done
shift $((OPTIND - 1))

# setup

# sanity
[[ ${ARGs} -gt 0 ]] && [[ $# -ne ${ARGs} ]] && abort "This script expects ${ARGs} arguments"
[[ -n ${RLIB} ]] && [[ ${BASE} -eq 1 ]] && abort "The options -b and -l are mutually exclusive"
[[ -n ${RLIB} ]] && [[ ! -d ${RLIB} ]] && abort "The directory ${RLIB} does not exist"
[[ -n ${RLIB} ]] && RLIB="-B ${RLIB}:/usr/local/lib/R/library"

# cmds container
cmds=()

# end of the boilerplate, logic goes below - instead of running cmds, add them to the cmds list
script=${1} && shift
cmds+=("apptainer exec -B /mnt:/mnt $RLIB $CONTAINER Rscript --vanilla $script --args $@
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
