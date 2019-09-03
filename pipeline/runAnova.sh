#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL

set -ex

## TODO add an option to make the name optional

usage(){
  echo >&2 \
  "
  This script expects three arguments: the transposed matrix, the gene list and the output dir.
  The transposed matrix should have no column nor row names.
  "
  exit 1
}

if [ $# != 3 ]; then
  echo "ERROR: This script expect three arguments"
  usage
fi

if [ ! -f $1 ]; then
  echo "ERROR: The first argument should be an existing file"
  usage
fi

if [ ! -f $2 ]; then
  echo "ERROR: The second argument should be an existing file"
  usage
fi

if [ ! -d $3 ]; then
  echo "ERROR: The third argument should be an existing directory"
  usage
fi

# create the structure
cd $3
ln -sf $1 NetworkE_expression_data.tsv
ln -sf $2 NetworkE_chip_features.tsv
nrow=`head -1 $1 | wc -w`

~bastian/Git/geneNetworkR/src/anova/anova NetworkE 1 $nrow
