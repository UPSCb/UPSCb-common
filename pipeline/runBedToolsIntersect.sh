#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## we get one dir and one file as input
export USAGETXT="Usage: $0 <a file> <b file> <out dir> [bed intersect option]"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# test the param
if [ "$#" -lt 3 ]; then
  abort "This function requires 3 arguments"
fi

if [ ! -f $1 ]; then
    abort "The first argument needs to be an existing file"
fi
a=$1
shift

if [ ! -f $1 ]; then
    abort "The second argument needs to be an existing file"
fi
b=$1;
shift;

if [ ! -d $1 ]; then
    abort "The third argument needs to be an existing directory"
fi
dir=$1;
shift;

# combine the filename for the output
outfile=$dir/`basename ${a%.*}`-`basename ${b%.*}`.tsv

## get the intersct results
bedtools intersect $@ -a $a -b $b > $outfile
