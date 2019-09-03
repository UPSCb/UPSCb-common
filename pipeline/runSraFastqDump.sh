#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mem=6GB
#SBATCH --mail-type=ALL

## stop on error and be verbose
set -ex

## usage
usage(){
  echo >&2 \
  "Usage: $(basename $0) <file to be converted> <output directory>"
  exit 1
}

## check number of arguments
if [ $# -lt 2 ]; then
   echo "This script takes two arguments"
   usage
fi

## check input file 
if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing sra file, please verify file path."
    usage
fi

## check output directory
if [ ! -d $2 ]; then
	echo "The second argument needs to be an existing output directory."
	usage
fi

## start
fastq-dump --gzip --split-3 -O $2 $1
# --split-3: create 3 files, one for forward reads, one for reverse, and one for unpaired reads.
# --I --split-files: produces two fastq files (--split-files) containing ".1" and ".2" read suffices (-I) for paired-end data
# --gzip: compress output using gzip
# - O: output directory
