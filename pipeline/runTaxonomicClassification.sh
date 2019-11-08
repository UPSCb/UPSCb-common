#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

# be safe (-e stop on error; -u stop if undefined variable, -x be verbose)
set -eux

# modules
module load bioinfo-tools vsearch/2.13.0

# load some functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# defaults
CUTOFF=0.9

# usage
USAGETXT=\
"
Usage: runVsearchTaxonClassification.sh <ITS2.fa> 

TODO: think that we might want a single file resulting from the dereplication, clustering, etc.
"

# process the arguments
if [ $# -ne 4 ]; then
         abort "This script expects 3 arguments"
fi 
     
if [ ! -f $1 ]; then
         abort "The first argument needs to be a file"
fi
     
if [ ! -f $2 ]; then
         abort "The second argument needs to be a file"
fi
     
if [ ! -d $3 ]; then
         abort "The third argument needs to be a directory"
fi

# run
vsearch --sintax $1 --db $2 --sintax_cutoff $CUTOFF --tabbedout $3/${1/.fa/.tsv} 

# think of compressing the output
