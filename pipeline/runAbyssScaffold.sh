#!/bin/bash

set -ex

## defaults
## set k to 25 for first run, increase to be more strict; k is the alignment length (KMER)
## n is the number of support for the edges, 1 is perfect for transcript; run incrementally for e.g. PacBio (EDGE)
PROC=32
KMER=25
EDGE=1

## usage
usage(){
    echo >&2 \
	"
        Usage: $0 [options] <BWA genome index> <transcripts fasta> <outdir>

        Options: -k kmer size; used for the bwa mem alignment; default to 25. Increase for stricter scaffolding.
                 -n edge support; default to 1 for transcripts; increase when using longer reads libraries instead
                 -p the number of threads to use in parallel
"
    exit 1
}

## get the options
while getopts k:n:p: option
do
    case "$option" in
	k) KMER=$OPTARG;;
	n) EDGE=$OPTARG;;
	p) PROC=$OPTARG;;
	\?) ## unknown flag
	    usage;;
    esac
done
shift `expr $OPTIND - 1`

## check the input
if [ $# != 3 ]; then
    echo "This script requires 3 arguments."
    usage
fi

if [ ! -f $1.bwt ]; then
    echo "The first argument needs to be the prefix - including the full path - of your BWA genome index"
    usage
fi
contigs=$1

if [ ! -f $2 ]; then
    echo "The second argument is the full path to your transcript/long reads fasta file "
    usage
fi
reads=$2

if [ ! -d $3 ]; then
    echo "The third argument is the full path to the desired output directory"
    usage
fi
out=$3

echo "Setting up"
if [ ! -d $out ]; then
    mkdir -p $out
fi

## get the basename
bnam=`basename ${reads//.f*a/}`
output=$out/$bnam

## 
echo "Aligning"
if [ ! -f $output.sam.gz ]; then
    bwa mem -a -t$PROC -S -P -k$KMER $contigs $reads | gzip > $output.sam.gz
fi

echo "Creating the graph"
if [ ! -f $output.dist.dot ]; then
    abyss-longseqdist -k$KMER $output.sam.gz | grep -v "l=" > $output.dist.dot
fi

echo "Scaffolding"
if [ ! -f $output.path ]; then
    abyss-scaffold -v -k$KMER -s200- -n$EDGE $contigs $output.dist.dot > $output.path
fi

echo "Consensus calling"
if [ ! -f ${output}2.path ]; then
    PathConsensus -v -k$KMER -p1 -s ${output}2.fa -g ${output}2.adj -o ${output}2.path $contigs $contigs $output.path
fi

echo "Merging"
cat $contigs ${output}2.fa | MergeContigs -v -k$KMER -o $output.fa - ${output}2.adj ${output}2.path

echo "Done"