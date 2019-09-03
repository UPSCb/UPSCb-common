#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL

## stop on error but be verbose
set -e
set -x

##
#  Run mmseq
##
## Usage:  sh runMmseq.sh fwd-fastq rv-fastq bowtie-index fasta-reference outputFolder

## executable
module load bioinfo-tools
module load mmseq

## arguments
if [ $# != 5 ]; then
   echo "This script takes four arguments: two fastq files, one bowtie index prefix and the output directory"
   exit 1
fi

## input files
if [ ! -f $1 ]; then
	echo "The first argument needs to be an existing fastq.gz file"
	exit 1
fi

if [ ! -f $2 ]; then
	echo "The second argument needs to be an existing fastq.gz file"
	exit 1
fi

## bowtie index
if [ ! -f $3.1.ebwt ]; then
	echo "The third argument needs to be an existing bowtie index (without the .1.ebwt extension)"
	exit 1
fi

## fasta ref
if [ ! -f $4 ]; then
	echo "The forth argument needs to be an existing fasta file"
	exit 1
fi

## output dir
if [ ! -d $5 ]; then
	echo "The forth argument needs to be an existing output directory." 
fi

## create the outfile name
outfile=$5/`basename ${1//_?[1,2]?.f*q.gz/}`

## IMPORTANT note
#--fullref (bowtie) and -m (bam2hits) are needed here as we do not have an EMSEMBL formatted header line
# see https://github.com/eturro/mmseq/#fastas-with-other-header-conventions
# TODO this would need adjustments to the script logic to make it parameterable!

## start
## run bowtie
bowtie -a --best --strata -S -m 100 -X 500 --chunkmbs 256 -p 16 $3 -1 <(gzip -dc $1) -2 <(gzip -dc $2) --fullref | samtools view -F 0xC -bS - | samtools sort -n - $outfile.bam

## then run bam2hits
bam2hits -m "(\S+)\s*.*" 1 1 $4 $outfile.bam > $outfile.hits

## then run
mmseq $oufile.hits $outfile

