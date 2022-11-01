#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -A facility
#SBATCH -J CHANGEME
#SBATCH --mail-usage=END,FAIL

# sanity
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Synopsis [options] $0 <ARGS>
  
  Options:
    -c change me

  Note: some notes
"

# vars
key="value"

# options
while getopts cC: option
do
    case "$option" in
    c) key="changeme-the-key-too";;
    C) key="changeme-the-key-too${OPTARG}";;
    \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# checks
[[ $# -ne CHANGEME ]] && abort "This script expects CHANGEME arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file path to a singularity container"
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"
[[ ! -f  ]] && abort "The CHANGEME argument needs to be an existing file path to CHANGEME"
[[ ! -d  ]] && abort "The CHANGEME argument needs to be an existing directory"

# run
singularity exec $1 \
CHANGEMETOOL CHANGEMECMDLINE
