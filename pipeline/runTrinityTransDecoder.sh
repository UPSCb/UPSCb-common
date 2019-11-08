#!/bin/bash -l
#SBATCH -p core
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00

# stop on error, be verbose and expand the commands
set -e -x

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
USAGETXT=\
"
	Usage: runTrinityTransDecoder.sh <Trinity.fasta> <out dir>
"

# Check
if [ $# -ne 2 ]; then
    echo "This function needs 2 arguments"
    usage
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the trinity fasta filepath"
fi

if [ ! -d $2 ]; then
    abort "The second argument (output dir) needs to be an existing directory"
fi

# run
cd $2

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
/usr/local/src/TransDecoder/TransDecoder.LongOrfs -t $1

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-trinotate-berlin2018.simg \
/usr/local/src/TransDecoder/TransDecoder.Predict -t $1
