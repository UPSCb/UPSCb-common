#!/bin/bash

# safeguards
set -ex

# project vars
account=SNIC2020-5-218
mail=nicolas.delhomme@umu.se

# directory
EXEC=/pfs/nobackup/home/b/bastian/seidr-devel/build/

# source
source functions.sh
source $EXEC/sourcefile

# Variables
chunkSize=350

# List of inference methods if we want to apply the chunk strategy to another
#inference=(aracne clr genie3 el-ensemble llr-ensemble narromi pcor pearson plsnet spearman svm-ensemble tigress)
inf=el-ensemble

# kk has 28 per node
# parameters
queueParams="-n 1 -c 28 -t 2-00:00:00"
commandParams="-B $chunkSize -O "'$SLURM_CPUS_PER_TASK'

# Small dataset params
smallDSet=0
minProp=80
maxProp=100

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

if [ $(wc -l $2 | cut -d" " -f1) -ne 1 ]; then
  abort "The second file should have only one line"
fi

# Setup
# create the output dir
mkdir -p results/$inf

# Read the gene list as an array
read -r -a GENEIDS <<< $(cat $2)
len=${#GENEIDS[@]}

# Extract the number of samples
nsamples=$(wc -l $1 | cut -d" " -f1)

# Extend the command line for small datasets
if [ $smallDSet -eq 1 ]; then
	smallDSetParam="-x $(expr $(expr $nsamples '*' $minProp) '/' 100) -X $(expr $(expr $nsamples '*' $maxProp) '/' 100)"
	commandParams="$commandParams $smallDSetParam"
fi

# Run
# Create chunks and iterate the submissions
for ((i=0;i<len;i+=$chunkSize)); do

  if [ ! -f results/$inf/$inf-$i.tsv ]; then
    printf "%s\n" "${GENEIDS[@]:$i:$chunkSize}" > gset-$i.txt

    srv=
    if [ -f results/$inf/$inf-$i.json ]; then
      srv="--resume-from results/$inf/$inf-$i.json"
    else
      srv="--save-resume results/$inf/$inf-$i.json"
    fi

    echo "#!/bin/bash" > results/$inf/$inf-$i.sh
    echo "unset OMP_NUM_THREADS" >> results/$inf/$inf-$i.sh
    echo "srun $EXEC/$inf $commandParams -i $1 -g $2 -t gset-$i.txt $srv -o results/$inf/$inf-$i.tsv" >> results/$inf/$inf-$i.sh
    sbatch --mail-type=ALL --mail-user=$mail -A $account -J $inf-$i -e results/$inf/$inf-$i.err -o results/$inf/$inf-$i.out $queueParams results/$inf/$inf-$i.sh
  fi
done

