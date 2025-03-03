#!/bin/bash -l
#SBATCH -p main -n 1
#SBATCH --time=1-00:00:00
#SBATCH --job-name="suppa2"
#SBATCH --mail-type=END,FAIL

# sanity
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Synopsis [options] $0 <suppa2 singularity container> <gtf file> <count tsv> <out dir>
"

# checks
[[ $# -ne 4 ]] && abort "This script expects four arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file path to a singularity container"
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"
[[ ! -f $2 ]] && abort "The second argument needs to be an existing file path to a gtf annotation file"
[[ ! -f $3 ]] && abort "The third argument needs to be an existing file path to a count tsv file"
[[ ! -d $4 ]] && abort "The fourth argument needs to be an existing directory"

# run
singularity exec $1 \
suppa.py psiPerIsoform \
-g $2 \
-o $4/$(basename ${3/.tsv/})_$(basename ${2/.gtf/}) \
-e $3
