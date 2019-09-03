#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --mail-type=ALL

## stop on error but be verbose
set -ex

##
#  Run fastQC 
usage(){
  echo >&2 \
  "Usage: $(basename $0) <outputFolder> <file> [file] ..."
  exit 1
}


## sanity checks
## executable
## are we on UPPMAX
if [ ! -z $SLURM_SUBMIT_DIR ]; then
	module load bioinfo-tools FastQC
##	echo "Running on UPPMAX"
else
##	echo "Running locally"
	fastqc=`which fastqc`
	if [ "$?" == "1" ]; then
		echo "please install fastqc before running this script or add it to your PATH"
		usage
	fi

	if [ ! -f $fastqc -a ! -x $fastqc ]; then
		echo "your fastQC does not appear to be an executable file"
		usage
	fi
fi

## arguments
if [ $# -lt 2 ]; then
   echo "This script takes two arguments"
   usage
fi

## input file
if [ ! -d $1 ]; then
	echo "The first argument needs to be an existing output directory."
	usage
fi
out=$1
shift

## output dir
for f in $@; do
  if [ ! -f $f ]; then
	echo "The second (and successive) arguments needs to be an existing fastq (optionally gz) file"
	usage
  fi
done

## start
fastqc --noextract --outdir $out $@ -t $SLURM_CPUS_ON_NODE
