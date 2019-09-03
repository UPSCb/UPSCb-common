#!/bin/bash -l
#SBATCH -p core
#SBATCH -c 1
#SBATCH --mail-type=ALL

set -ex

# usage
usage(){
echo >&2 \
"
	Usage: $0 <expression matrix> <metadata table> <output dir>
	
	The expression matrix and metadata table have to be in tab delimited format
	and may be gzipped.
"
	exit 1
}

# check
if [ $# != 3 ]; then
  echo "This script expects 3 arguments"
  usage
fi

if [ ! -f $1 ]; then
  echo "The first argument must be a file"
  usage
fi

if [ ! -f $2 ]; then
  echo "The second argument must be a file"
  usage
fi

if [ ! -d $3 ]; then
  echo "The third argument must be a directory"
  usage
fi

# get the exec
module load R
exeR=`Rscript -e 'cat(system.file("R","geneNetworkR-preparation.R",package="geneNetworkR"))'`


# run with knitr
Rscript -e "library(knitr); spin('$exeR')" -e $1 -m $2 -f $3

