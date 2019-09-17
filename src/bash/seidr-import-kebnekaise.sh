#!/bin/bash

# safeguards
set -ex

# project vars
account=SNIC2018-3-61
mail=nicolas.delhomme@umu.se

# source
source functions.sh

# modules
#module load bioinfo-tools seidr-devel
#export PATH=/pfs/nobackup/home/b/bastian/seidr/build:$PATH
source /pfs/nobackup/home/b/bastian/seidr/build/sourcefile

# Variables

inference=(aracne clr elnet genie3 llr-ensemble mi narromi pcor pearson plsnet spearman tigress)

CPUs=14

# usage
USAGETXT=\
"
$0 <inference.tsv> <genes.tsv>

Methods to be imported: ${inference[@]}

"

# input: gene and data
if [ $# -ne 2 ]; then
  abort "This script expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the inference result file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}
  if [ ! -f results/$inf/$inf.sf ]; then
    
    /pfs/nobackup/home/b/bastian/seidr/build/scripts/generate_import_script.py \
    -i results/$inf/$inf.tsv -g genes.tsv -c $CPUs > results/$inf/${inf}-import.sh
    
    # echo "#!/bin/bash" > results/$inf/$inf.sh
    # echo "unset OMP_NUM_THREADS" >> results/$inf/$inf.sh
    # echo "srun ${command[$i]} -i $1 -g $2 -o results/$inf/$inf.tsv" >> results/$inf/$inf.sh
    sbatch --mail-type=ALL --mail-user=$mail -A $account -J import-$inf \
    -e results/$inf/${inf}-import.err -o results/$inf/${inf}-import.out results/$inf/${inf}-import.sh
  
    # ${arguments[$i]} 
  fi
done

