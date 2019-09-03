#!/bin/bash

# safeguards
set -ex

# project vars
account=u2018015
mail=nicolas.delhomme@umu.se

# check
if [ -z $UPSCb ]; then
  echo "The UPSCb environment variable needs to be set to your UPSCb Git checkout path"
  exit 1
fi

# source
source $UPSCb/src/bash/functions.sh

# modules
module load bioinfo-tools seidr-devel
#export PATH=/pfs/nobackup/home/b/bastian/seidr/build:$PATH
#source /pfs/nobackup/home/b/bastian/seidr/build/sourcefile

# Variables

inference=(aracne clr elnet genie3 llr-ensemble narromi pcor pearson plsnet spearman tigress)

# 14 workers on 2 nodes (kk has 28 per node)
# narromi is not thread safe, hence -c 1
#default="-n 2 -c 14 -t 1-00:00:00"
#arguments=([0]=$default [1]=$default [2]="-n 2 -c 28 -t 2-00:00:00" [3]=$default [4]="-n 28 -c 1 -t 1-00:00:00" [5]="-n 7 -t 12:00:00" [6]="-n 7 -t 12:00:00" [7]=$default [8]="-n 7 -t 12:00:00" [9]="-n 3 -c 28 -t 3-00:00:00")
#parallel="-O "'$SLURM_CPUS_PER_TASK'
#command=([0]="mi -m aracne "$parallel [1]="mi -m CLR "$parallel [2]="genie3 "$parallel [3]="llr-ensemble "$parallel [4]="narromi "$parallel [5]="pcor" [6]="correlation -m pearson" [7]="plsnet "$parallel [8]="correlation -m spearman" [9]="tigress "$parallel)
CPUs=14

# usage
USAGETXT=\
"
$0 <genes.tsv>

Methods to be imported: ${inference[@]}

"

# input: gene and data
if [ $# -ne 1 ]; then
  abort "This script expects 1 argument"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}

  if [ -f results/$inf/$inf.tsv ]; then

    $UPSCb/src/python/generate_import_script.py \
    -i results/$inf/$inf.tsv -g genes.tsv -c $CPUs  > results/$inf/${inf}-import.sh

    sbatch --mail-type=ALL --mail-user=$mail -A $account -J import-$inf -e results/$inf/${inf}-import.err -o results/$inf/${inf}-import.out results/$inf/${inf}-import.sh
  fi

done

