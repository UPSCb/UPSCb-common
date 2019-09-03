#!/bin/bash -l
#SBATCH -p core
#SBATCH -c 1
#SBATCH --mail-type=ALL

set -ex

# usage
usage(){
echo >&2 \
"
	Usage: $0 [options] <output dir> <tool>

	Note: The tool is one of the following:
	 Anova
	 CLR
	 GeneNet
	 GENIE3
	 Narromi
	 Pearson
	 Spearman
	 TIGLM
	 TIGRESS

	 Options:
                -c the number of cores to parallelise over
"
	exit 1
}

# define global vars
CPU=1

# manage options
while getopts c: option
do
  case "$option" in
	    c) CPU=$OPTARG;;
	    \?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

# check
if [ $# != 2 ]; then
  echo "This script expects 2 arguments"
  usage
fi

if [ ! -d $1 ]; then
  echo "The first argument must be a directory"
  usage
fi

if [ ! -f $1/Data/sexp.rda ]; then
  echo "The directory has to be a valid geneNetworkR directory."
  echo "Run geneNetworkRPreparation.R first."
  usage
fi

# We do not check $2 (the tool) as this is done in the Rscript

# get the exec
module load R
exeR=`Rscript -e 'cat(system.file("R","geneNetworkR-run.R",package="geneNetworkR"))'`

# run with knitr
Rscript -e "library(knitr); spin('$exeR')" -c $CPU -f $1 -t $2

