#!/bin/bash
#SBATCH -p core -n 1
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
# stop on error
set -ex

# usage
export USAGETXT=\
"
Usage: $0 <FASTA> <OUTFILE>
"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# checks
if [ $@ -ne 2 ]; then
  abort "This script expects 2 arguments."
fi

if [ ! -f $1 ]; then
  abort "The fasta file does not exist"
fi

if [ ! -d $(dirname $outfile) ]; then
  abort "The directory of the output file does not exist"
fi

isExec muscle

# run
muscle -in $1 -out $2 -diags -maxiters 1
