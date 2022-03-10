#!/bin/bash
#SBATCH -p core
#SBATCH -n 2
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL
#SBATCH -t 12:00:00

## stop on error but be verbose
set -eux

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: $(basename $0) <singularity container> <analysis directory> <output directory>
"

## arguments
[[ $# -ne 3 ]] && abort "This script takes three arguments"


[[ ! -f $1 ]] && abort "The first argument needs to be an existing mutiqc singularity container file."

[[ ! -d $2 ]] && abort "The second argument needs to be an existing analysis directory."

[[ ! -d $3 ]] && abort "The third argument needs to be an existing output directory."

## start
singularity exec --bind /mnt:/mnt $1 multiqc -o $3 $2
