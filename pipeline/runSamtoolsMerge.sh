#!/bin/bash -l
#SBATCH -p main
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL

# be verbose and print
set -ex

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# test
#isEnvVarSet $UPSCb

# usage
USAGETXT=\
"
Usage: runSamtoolsMerge.sh <samtools singularity container> <output bam file> <input bamfile 1> <input bamfile 2> ... <input bamfile n>
"

[[ $# -lt 4 ]] && abort "The script expects at least four arguments"

[[ ! -f $1 ]] && abort "The singularity container needs to be a file"
singularity=$1
shift

[[ ! -d $(dirname $1) ]] && abort "The output directory of the output file does not exist."
out=$1
shift

for f in $@; do
  [[ ! -f $f ]] && abort "The input BAM $f does not exist"
done

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

#run samtools merge
singularity exec $singularity samtools merge $out $@
