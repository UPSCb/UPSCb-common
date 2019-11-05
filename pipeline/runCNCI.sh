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
	Usage: runCNCI.sh <Trinity.fasta> <out dir>
	
	Options:
	            -f    input file
	            -o    output file 
	            -p    (parallel) assign the running CUP numbers
	            -m    (model) assign the classification model ("re" for vertebrates; "pl" for plants)
	            
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

# run CNCI
cd $2

singularity exec --bind /mnt:/mnt \
/mnt/picea/projects/singularity/delhomme-upscb-lncrna.simg CNCI.py \
-p 12 -m pl -f $1 -o $2


