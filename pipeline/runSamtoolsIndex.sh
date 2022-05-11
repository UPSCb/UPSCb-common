#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=END,FAIL

# stop on error
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# Vars
CSI=

# a usage function
USAGETXT=\
"
Usage: $0 [options] <samtools singularity container> <bam file>

Options: -c generate CSI index (for large, e.g. the spruce, genomes)

" 

# get the options
while getopts ch option
do
        case "$option" in
        c) CSI="-c";;
	      h) usage;;
		    \?) ## unknown flag
		    abort "unknown option";;
        esac
done
shift `expr $OPTIND - 1`

# safeguards
[[ $# != 2 ]] && abort "This function takes two arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be the singularity container file"

[[ ! -f $2 ]] && abort "The first argument needs to be an existing bam file"

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# create the index
singularity exec $1 samtools index $CSI $2
