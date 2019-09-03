#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e

##
echo Checking

## we get two dir and two files as input
if [ $# != 2 ]; then
    echo "This function takes two files as arguments"
    echo "Usage: sbatch runRePair.sh <forward fastq> <reverse fastq>"
    echo "Note: The UPSCb env. var. needs to be set to your UPSCb Git checkout dir."
    exit 1
fi

if [ ! -f $1 ]; then
    echo "The second argument needs to be an existing fastq file"
fi

if [ ! -f $2 ]; then
    echo "The third argument needs to be an existing fastq file"
fi

##
echo Pairing

## check the file order
## SLURM_SUBMIT_DIR is the directory from which sbatch was invoked.
## in our case it should always be in the project pipeline dir, so 3 dirs
## down the hierarchy to a common parent directory
Rscript $UPSCb/src/R/rePairFastq.R -f $1 -r $2 -v

## 
echo Gzipping

## compress the output files
printf "%s\0%s" ${1/1./_paired_1.} ${2/2./_paired_2.} ${1/1./_single_1.} ${2/2./_single_2.} | xargs -0 -I {} -P 4 gzip -f {}

##
echo Done
