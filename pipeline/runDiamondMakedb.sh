#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=FAIL
#SBATCH -t 1:00:00

set -eux

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
  Usage: $0 <singularity diamond container><fasta> <indexName>
"

[[ $# -ne 3 ]] && abort "This script expects 3 arguments"

[[ ! -f $1 ]] && abort "The singularity container file does not exist"

[[ ! -f $2 ]] && abort "The input fasta file does not exist"

[[ ! -d $(dirname $3) ]] && abort "The output directory for the index does not exist"

[[ ${SINGULARITY_BINDPATH:-1} -eq 1 ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

singularity exec $1 diamond makedb --in $2 -d $3
