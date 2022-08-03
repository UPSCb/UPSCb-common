#!/bin/bash -l
#SBATCH --mail-type=END,FAIL
#SBATCH -p core -w picea
#SBATCH --mem=264GB
#SBATCH -t 12:00:00

# failsafe
set -eu

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
    runBwaIndex.sh <bwa singularity container> <Genome Fasta> <Output Dir>
"

# sanity
[[ $# -ne 3 ]] && abort "The script expects three arguments."
[[ ! -f $1 ]] && abort "BWA singularity container not found"
[[ ! -f $2 ]] && abort "FASTA input not found"
[[ ! -d $3 ]] && abort "OUTPUT directory not found"

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# prep
ln -sf $2 $3
BNAM=$(basename $2)
singularity exec $1 bwa index $3/$BNAM
