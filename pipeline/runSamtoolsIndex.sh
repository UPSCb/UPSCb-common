#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL

## stop on error
set -e

## modules
module load bioinfo-tools samtools

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 <fasta file>
" 
    exit 1
}

## we get one file as input
if [ $# != 1 ]; then
    echo "This function takes one fasta file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing fasta file"
    usage
fi

## create the index
samtools faidx $1

