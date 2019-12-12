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
	Usage: runPLEK.sh <Trinity.fasta> <out dir>
	
	Options:
	            -thread       number of threads for running the PLEK program
	            -minlength    the minimum length of sequences
	            -isoutmsg     output messages to screen or not
	            -isrtempfile  remove temporary files or not
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

# run PLEK
cd $2

# output filename
fnam=$(basename ${1%.*})
singularity exec --bind /mnt:/mnt \
/mnt/picea/projects/singularity/delhomme-upscb-lncrna.simg PLEK.py \
-minlength 200 -thread 12 -p 12 -isoutmsg 1 -fasta $1 -out $2/$fnam.txt


