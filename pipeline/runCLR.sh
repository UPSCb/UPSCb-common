#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL

set -ex

usage(){
  echo >&2 \
  "
  This script expects two arguments: the transposed matrix and the output dir.
  The transposed matrix should have no column nor row names.
  "
  exit 1
}

if [ $# != 2 ]; then
  echo "ERROR: This script expects two arguments"
  usage
fi

if [ ! -f $1 ]; then
  echo "ERROR: The first argument should be an existing file"
  usage
fi

if [ ! -d $2 ]; then
  echo "ERROR: The second argument should be an existing directory"
  usage
fi

# global vars
gp=$UPSCb/src/c/genepair
rs=$UPSCb/src/NetworkCrowd/CLR/CLR.R

# load module
module load R

# create the structure
cd $2
ln -sf $gp .
ln -sf $rs .

# run
Rscript --vanilla $rs --data $1
