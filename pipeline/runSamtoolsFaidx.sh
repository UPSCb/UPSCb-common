#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL

## stop on error
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## a usage function
USAGETXT=\
"
Usage: $0 <samtools singularity container> <fasta file>
"

# safeguards
[[ $# != 2 ]] && abort "This function takes two arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be the singularity container file"

[[ ! -f $2 ]] && abort "The first argument needs to be an existing bam file"

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

## create the index
singularity exec $1 samtools faidx $2
