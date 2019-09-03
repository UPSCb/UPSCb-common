#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL

## stop on error but be verbose
set -ex

usage(){
echo >&2 \
"
	Usage: $0 [option] <out dir> <genome fasta> <in bam>
"
	exit 1
}

## executable
module load bioinfo-tools
module load mmseq

## arguments
if [ $# != 3 ]; then
   echo "This script takes three arguments: one out dir, one reference fasta file and one input bam file"
   exit 1
fi

## input files
if [ ! -d $1 ]; then
	echo "The first argument needs to be an existing directory"
	usage
fi

if [ ! -f $2 ]; then
	echo "The second argument needs to be an existing fastq(.gz) file"
	usage
fi

## bowtie index
if [ ! -f $3 ]; then
	echo "The third argument needs to be an existing bam file"
	usage
fi

## create the outfile name
outfile=$1/`basename ${3//.bam/}`

## start
## then run bam2hits
bam2hits -m "(\S+)\s*.*" 1 1 $2 $3 > $outfile.hits

## then run
mmseq $outfile.hits $outfile

