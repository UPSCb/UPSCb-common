#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## load the modules
module load bioinfo-tools
module load BEDTools
module load samtools

## we get one dir and one file as input
usage(){
    echo >&2 \
    "Usage: $0 <bam file> <fastq fwd file> <fastq rev file>"
    exit 1
}

if [ $# != 3 ]; then
  echo "This function requires 3 arguments"
  usage;
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing file"
    usage;
fi

## samtools
samtools sort -n $1 ${1//.bam/.nsorted}

## extract the fastq files
bedtools bamtofastq -i ${1//.bam/.nsorted}.bam -fq $2 -fq2 $3

## cleanup
rm ${1//.bam/.nsorted}.bam
gzip -f $2
gzip -f $3
