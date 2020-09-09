#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem=4GB
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage 
USAGETXT=\
"
	Usage: $0 <in.bam>
"

## we get one file as input
if [ $# != 1 ]; then
    echo "This function takes one file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing bam file"
fi

## define the output file
new=${1//.bam/.stats}

## get the coverage table
samtools flagstat $1 > $new
