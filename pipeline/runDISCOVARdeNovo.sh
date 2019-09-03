#!/bin/bash -l
#SBATCH --mail-type=ALL
module load jemalloc gcc
usage () 
{
	echo "runDISCOVARdeNovo.sh <fastq1,fastq2,fastq3,fastq4...> <out_dir> <Discovar options>"
	echo
	echo "Note: This script expects the input FASTQ (or SAM/BAM) reads in a"
	echo "comma separated list just as DISCOVAR does"
	exit 1
}

#Argn check

if [ $# -lt 2 ]; then
	usage
fi

# File and dir errors
for f in $(echo $1 | tr "," " "); do
	if [ ! -f $f ];
		then echo "File $f not found. Exiting"
		usage
	fi
done

if [ ! -d $2 ]; then
	echo "Directory $2 note found. Exiting"
	usage
fi

# Arg assignments
INF=$1
OUTD=$2

shift 2

EXTRA_ARGS=$@


DiscovarDeNovo READS=$INF OUT_DIR=$OUTD $EXTRA_ARGS
