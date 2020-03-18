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
chunkSize=1000

# List of inference methods if we want to apply the chunk strategy to another
#inference=(aracne clr genie3 el-ensemble llr-ensemble narromi pcor pearson plsnet spearman svm-ensemble tigress)
inf=el-ensemble

# kk has 28 per node
# parameters
queueParams="-n 2 -c 14 -t 2-00:00:00"
commandParams="-O "'$SLURM_CPUS_PER_TASK'

# usage
USAGETXT=\
"
$0 <expression-matrix.tsv> <genes.tsv> 
"

# Validation
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

# Setup
# create the output dir
mkdir -p results/$inf

# Read the gene list as an array
read -r -a GENEIDS <<< $(cat $2)
len=${#GENEIDS[@]}

# Run
# Create chunks and iterate the submissions
for ((i=0;i<len;i+=$chunkSize)); do
  
  if [ ! -f results/$inf/$inf-$i.tsv ]; then
    printf "%s\n" "${GENEIDS[@]:$i:$chunkSize}" > gset-$i.txt
  
    echo "#!/bin/bash" > results/$inf/$inf-$i.sh
    echo "unset OMP_NUM_THREADS" >> results/$inf/$inf-$i.sh
    echo "srun $inf $commandParams -i $1 -g $2 -t gset-$i.txt -o results/$inf/$inf-$i.tsv" >> results/$inf/$inf-$i.sh
    sbatch --mail-type=ALL --mail-user=$mail -A $account -J $inf-$i -e results/$inf/$inf-$i.err -o results/$inf/$inf-$i.out $queueParams results/$inf/$inf-$i.sh
  fi
done

