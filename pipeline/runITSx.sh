#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

# be safe (-e stop on error; -u stop if undefined variable, -x be verbose)
set -eux

# load some functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

CPU=12

USAGETXT=\
"
Usage: runITSx.sh <mergedpairs fasta file> <forward fasta file> <output directory> <ITSx HMMs directory)
"

# process the arguments
if [ $# -ne 4 ]; then
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
     
if [ ! -d $4 ]; then
     abort "The forth argument needs to be a directory"
fi
     
# merge the input
input=$3/$(basename ${1/.fa/-tmp.fa})
cat $1 $2 > $input
     
# run
ITSx -i $input -o ${input/-tmp.fa/-ITSx.fa} -p $4 --cpu $CPU
     
# clean
rm $input
