#! /bin/bash

## stop on error
set -e

## be verbose and extend the commands
set -x

## usage
usage(){
echo >&2 \
"
	Usage: runBamSubset.sh <contig list file> <bam file>

	The <contig list file> should contain one contig name per line
"
	exit 1
}

if [ $# != 2 ]; then
    echo "This function takes two files as arguments"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing file listing the contigs"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing bam file"
    usage
fi

## get the out name
out=${2//.bam/_Subset}

## the awk one liner
samtools view -h $2 | awk '{if(NF==1){ctg[$1]++} else {if((x=index($1,"@")) > 0){if($1 == "@SQ"){sn=$2;sub("SN:","",sn);if (ctg[sn]==1){print $0}}else{print $0}} else {if (ctg[$3]==1){print $0}}}}' $1 - | samtools view -bS - | samtools sort - $out