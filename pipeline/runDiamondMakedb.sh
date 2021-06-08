#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=FAIL
#SBATCH -t 1:00:00

set -eux

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
  Usage: $0 <fasta> <indexDir>
"

[[ $# -ne 2 ]] && abort "This script expects 2 arguments"

[[ ! -f $1 ]] && abort "The input fasta file does not exist"

[[ ! -d $2 ]] && abort "The output directory does not exist"

diamond makedb --in $1 -d $2
