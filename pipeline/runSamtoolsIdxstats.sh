#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -e

## modules
module load bioinfo-tools
module load samtools/0.1.19

## we get one dir and one file as input
if [ $# != 2 ]; then
    echo "This function takes one directories and one bam file as arguments"
    echo "Usage: sbatch runSamtoolsIdxstats.sh <out dir> <in.bam>"
    exit 1
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing bam file"
fi
nam=`basename ${2//.bam/}`

## get the coverage table
samtools idxstats $2 | cut -f 1,3 > $1/$nam.txt

