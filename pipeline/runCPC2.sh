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
	Usage: runCPC2.sh <Trinity.fasta> <out dir>
	
	Options:
	            -i    input file
	            -o    output file 
	            -r    also check the reverse strand [Default: FALSE]
	            
"

# Check
if [ $# -ne 2 ]; then
    echo "This function needs 2 arguments"
    usage
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the trinity fasta filepath"
fi

if [ ! -d $(dirname $2) ]; then
    abort "The second argument (output dir) parent directory does not exist"
fi

# run CPC2
cd $(dirname $2)
docker run --rm -v /mnt:/mnt \
delhomme/upscb-lncrna CPC2.py \
-r TRUE -i $1 -o $2


