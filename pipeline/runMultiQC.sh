#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00

## stop on error but be verbose
set -ex

##
#  Run MultiQC 
usage(){
  echo >&2 \
  "Usage: $(basename $0) <analysis directory> <output directory>"
  exit 1
}

## arguments
if [ $# -lt 2 ]; then
   echo "This script takes two arguments"
   usage
fi

## input file
if [ ! -d $1 ]; then
	echo "The first argument needs to be an existing analysis directory."
	usage
fi

if [ ! -d $2 ]; then
	echo "The second argument needs to be an existing output directory."
	usage
fi

## start
multiqc -o $2 $1
