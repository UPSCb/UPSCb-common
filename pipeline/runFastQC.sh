#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --mail-type=END,FAIL

# fail on ERROR
set -eux

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# vars
OPTIONS="--noextract"
CPU=1

# usage
USAGETXT=\
"
 $0 <singularity image> <outputFolder> <fastq file>
"

## arguments
[[ $# -ne 3 ]] && abort "This script takes two arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity fastqc container file"
	
[[ ! -d $2 ]] && abort "The second argument needs to be an existing directory"

[[ ! -f $3 ]] && abort "The third argument needs to be an fastq file"

## start
singularity exec $1 fastqc --outdir $2 -t $CPU $OPTIONS $3
