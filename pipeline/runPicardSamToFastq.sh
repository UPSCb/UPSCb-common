#!/bin/bash -l
#SBATCH -t 6:00:00
#SBATCH -p main -n 2
#SBATCH --mail-type=FAIL

set -eu

# setup
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage  $0 <bam> <out>
"

# validity
[[ $# -ne 2 ]] && abort "This script expects 2 arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be a bam file"
[[ ! -d $2 ]] && abort "The second argument needs to be a directory"

# extract sample name
fnam=$(basename ${1/.bam/})

# extract fastq
java -jar $PICARD_ROOT/picard.jar SamToFastq \
-I $1 -F $2/${fnam}_1.fq -F2 $2/${fnam}_2.fq -FU $2/${fnam}.fq \
--VALIDATION_STRINGENCY LENIENT

# compress
gzip -f $2/${fnam}_1.fq
gzip -f $2/${fnam}_2.fq

# clean
rm $2/${fnam}.fq
