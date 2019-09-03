#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools star-fusion samtools

# OPTIONS
CPU=8

# usage function
usage(){
echo >&2 \
"
	Usage: $0 <star-fusion index> <star chimeric file> <out dir>
"
	exit 1
}

# check the arguments
if [ ! -d $1 ]; then
	echo "The star-fusion index dir: $1 does not exist"
	usage
fi

if [ ! -f $1/_ref_cdna.fasta.bowtie_idx.ok ]; then
	echo "The star-fusion index directory: $1 does not seem to be a valid index directory"
	usage
fi

if [ ! -f $2 ]; then
	echo "The star chimeric file: $2 does not exist"
	usage
fi

if [ ! -d $3 ]; then
  echo "The output directory: $3 does not exist"
  usage
fi

# run the commands
STAR-Fusion --genome_lib_dir $1 \
             -J $2 \
             --output_dir $3
