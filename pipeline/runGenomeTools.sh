#!/bin/bash
#SBATCH -p core
#SBATCH -n 2
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL
#SBATCH -t 12:00:00

# stop on error and undefined vars
set -eu

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: <singularity genometools> [options] <output gff file> <input gff file>
"
## arguments
[[ $# -ne 3 ]] && abort "This script takes three arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing genome tools singularity container file."
[[ ! -d $2 ]] && abort "The second argument needs to be the output gff filename."
[[ ! -d $3 ]] && abort "The third argument needs to be an existing input gff file."

## start
singularity exec $1 gt gff3 -force -tidy yes -addintrons yes -addids yes \
-fixregionboundaries yes -retainids yes -sort yes -checkids yes -o $2 $3
