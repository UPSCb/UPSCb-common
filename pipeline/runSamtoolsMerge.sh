#!/bin/bash -l

#load module(s)
module load bioinfo-tools samtools

# be verbose and print
set -ex

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# test
#isEnvVarSet $UPSCb

# usage
USAGETXT=\
"
Usage: runSamtoolsMerge.sh <output bam file> <input bamfile 1> <input bamfile 2>
"

# test if $2 - $3 are -f

if [ ! -d $(dirname $1)]; then
    abort "The output directory does not exist."
fi

out=$1
shift

#run samtools merge
samtools merge $out $@

samtools index $out
