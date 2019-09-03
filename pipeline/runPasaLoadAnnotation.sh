#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools pasa

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 <config file> <genome fasta file> <gff3 file>
" 
    exit 1
}

## we get three files as input
if [ $# != 3 ]; then
    echo "This function takes one config file, one fasta file and one gff3 file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing config file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing fasta file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing gff3 file"
    usage
fi

## execute
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi -c $1 -g $2 -P $3
