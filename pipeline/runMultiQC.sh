#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00

# stop on error but be verbose
set -eux

# load helpers
source functions.sh

#  Run MultiQC
USAGETXT=\
"Usage: $(basename $0) <singularity image> <analysis directory> <output directory>"

## arguments
[[ $# -lt 3 ]] &&  abort "This script takes three arguments"

## input file
[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity multiqc container file"

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# directory
[[ ! -d $2 ]] && abort "The second argument needs to be an existing analysis directory."

# output
[[ ! -d $3 ]] && abort "The third argument needs to be an existing output directory."

## start
singularity exec $1 multiqc -o $3 $2
