#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## load the modules
module load bioinfo-tools
module load BEDTools

## we get one dir and one file as input
usage(){
    echo >&2 \
    " Usage: $0 <a file> <b file> <out dir> [bed subtract option] 
    "
    exit 1
}

if [ "$#" -lt 3 ]; then
  echo "This function requires 3 arguments"
  usage;
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing file"    
    usage;
fi
a=$1
shift

if [ ! -f $1 ]; then
    echo "The second argument needs to be an existing file"
    usage;
fi
b=$1;
shift;

if [ ! -d $1 ]; then
    echo "The third argument needs to be an existing directory"
    usage;
fi
dir=$1;
shift;

# combine the filename for the output
outfile=$dir/`basename ${a%.*}`-`basename ${b%.*}`."${a##*.}"

## get the subtracted results
bedtools subtract $@ -a $a -b $b > $outfile
