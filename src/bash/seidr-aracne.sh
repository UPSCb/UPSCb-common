#!/bin/bash

# safeguards
set -ex

# project vars
account=facility
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

inf=aracne

# 14 workers on 2 nodes (kk has 28 per node)
# narromi is not thread safe, hence -c 1
arguments="-n 2 -c 14 -t 1-00:00:00"
command="mi -m ARACNE -O "'$SLURM_CPUS_PER_TASK'

# usage
USAGETXT=\
"
$0 <expression-matrix.tsv> <genes.tsv>
"

# input: gene and data
if [ $# -ne 2 ]; then
  abort "This script expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the expression matrix tab delimited file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
if [ ! -f results/$inf/$inf.tsv ]; then
  mkdir -p results/$inf
  echo "#!/bin/bash" > results/$inf/$inf.sh
  echo "unset OMP_NUM_THREADS" >> results/$inf/$inf.sh
  echo "srun $command -i $1 -g $2 -o results/$inf/$inf.tsv" >> results/$inf/$inf.sh
  sbatch --mail-type=ALL --mail-user=$mail -A $account -J $inf \
  -e results/$inf/$inf.err -o results/$inf/$inf.out $arguments results/$inf/$inf.sh
fi


