#!/bin/bash
#SBATCH -p all
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -t 2-00:00:00

# modules
module load R

# check
if [ -z $UPSCb ]; then
  "You need to have the UPSCb env. var. set to your UPSCb Git checkout dir"
fi

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage: $0 <R script to knit>
"

# sanity: this script expects one argument, the file to knit
if [ $# -ne 1 ]; then
  abort "This script expects one argument, the file to knit"
fi

# knit
Rscript -e "require(methods);rmarkdown::render(commandArgs(trailingOnly=TRUE)[1])" $1
