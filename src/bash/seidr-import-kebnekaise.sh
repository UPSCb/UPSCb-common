#!/bin/bash

# safeguards
set -ex

# project vars
account=SNIC2019-3-207
mail=nicolas.delhomme@slu.se

# source
source functions.sh

# modules
source /pfs/nobackup/home/b/bastian/seidr/build/sourcefile

# Variables
inference=(aracne clr elnet genie3 llr-ensemble mi narromi pcor pearson plsnet spearman tigress)

# additional parameters (elnet is done iteratively, so the format is not the expected one: a matrix, rather an edge list
arguments=([2]="-f el" [4]="-o llr-ensemble.sf")
CPUs=28
Time=1-00:00:00

# usage
USAGETXT=\
"
$0 <genes.tsv>

Methods to be import from: ${inference[@]}

"

# input: gene and data
if [ $# -ne 1 ]; then
  abort "This script expects 1 argument"
fi

if [ ! -f $1 ]; then
  abort "The argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}
  if [ ! -f results/$inf/$inf.sf ]; then

    if [ -f results/$inf/$inf.tsv ]; then
      /pfs/nobackup/home/d/delhomme/seidr/scripts/generate_import_script.py \
      -i results/$inf/$inf.tsv -g $1 -c $CPUs ${arguments[$i]} > results/$inf/${inf}-import.sh

      sbatch -t $Time --mail-type=ALL --mail-user=$mail -A $account -J import-$inf \
	-e results/$inf/${inf}-import.err -o results/$inf/${inf}-import.out results/$inf/${inf}-import.sh
    else
      echo "There is no tsv file for $inf"
    fi

    # ${arguments[$i]}
  fi
done


