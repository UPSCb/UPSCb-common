#!/bin/bash -l

## report error
set -e

## be verbose 
set -x

## usage
usage(){
echo >&2 \
"
     Usage: runBESST.sh <genome fasta> <alignment bam> <out dir>
"
exit 1
}

## check params
if [ $# != 3 ]; then
    echo "This script needs three parameters"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first parameter has to be the genome (scaffold) fasta file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second parameter should be a bam file"
    usage
fi

if [ ! -f $2.bai ]; then
    echo "The bam file should be indexed, i.e. a $1.bai file should exist."
    usage
fi

if [ ! -d $3 ]; then
    echo "The third argument has to be an existing output directory"
    usage
fi

# load modules
module load bioinfo-tools BESST

## run
runBESST -c $1 -f $2 -o $3
