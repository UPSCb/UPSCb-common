#!/bin/bash -l
#SBATCH -p main
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL

# be safe (-e stop on error; -u stop if undefined variable, -x be verbose)
set -eux

# load some functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: runVsearchMergePair.sh <forward fastq file> <reverse fastq file> <output directory>
"

# process the arguments
if [ $# -ne 3 ]; then
    abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be a file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be a file"
fi

if [ ! -d $3 ]; then
  abort "The third argument needs to be a directory"
fi

# run
fnam=$(basename ${1/_1.fastq.gz/})
vsearch --fastq_mergepairs $1 --reverse $2 --fastq_allowmergestagger \
--fastaout $3/$fnam.fa --fastaout_notmerged_fwd $3/${fnam}_1.fa --fastaout_notmerged_rev $3/${fnam}_2.fa

# think of compressing the output
