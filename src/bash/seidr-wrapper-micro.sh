#!/bin/bash

# safeguards
set -ex

# project vars
account=facility
mail=

# source
source functions.sh

# modules
module load bioinfo-tools clp seidr-devel

# Variables
resultDir=results

nodes=1

genie3Multiplier=2

tigressMultiplier=3

partition="rbx"
#partition="mpi"

correlationNCPUs=24

inference=(aracne clr genie3 llr-ensemble narromi pcor pearson plsnet spearman tigress tomsimilarity)

# which inference to run
run=(
  [0]=1
  [1]=1
  [2]=1
  [3]=1
  [4]=1
  [5]=1
  [6]=1
  [7]=1
  [8]=1
  [9]=1
  [10]=1)

# queue
default="-n $nodes -c $correlationNCPUs -t 2-00:00:00"

arguments=(
  [0]=$default
  [1]=$default
  #[2]="-n $(expr $nodes "*" genie3Multiplier) -c $(expr $correlationNCPUs "*" genie3Multiplier) -t 4-00:00:00"
  [2]="-n $nodes -c $(expr $correlationNCPUs "*" genie3Multiplier) -t 4-00:00:00"
  [3]=$default
  [4]=$default
  [5]=$default
  [6]=$default
  [7]=$default
  [8]=$default
  #[9]="-n $(expr $nodes "*" tigressMultiplier) -c $(expr $correlationNCPUs "*" tigressMultiplier) -t 4-00:00:00"
  [9]="-n $nodes -c $(expr $correlationNCPUs "*" tigressMultiplier) -t 4-00:00:00"
  [10]=$default)
  
#parallel="-O "'$SLURM_CPUS_PER_TASK'
parallel="-O $correlationNCPUs"
parallel2="-O "$(expr $correlationNCPUs "*" genie3Multiplier)
parallel3="-O "$(expr $correlationNCPUs "*" tigressMultiplier)
command=(
  [0]="mi -m ARACNE -M $resultDir/mi/mi.tsv "$parallel
  [1]="mi -m CLR -M $resultDir/mi/mi.tsv "$parallel
  [2]="genie3 "$parallel2
  [3]="llr-ensemble "$parallel
  [4]="narromi -m interior-point "$parallel
  [5]="pcor"
  [6]="correlation -m pearson"
  [7]="plsnet "$parallel
  [8]="correlation -m spearman"
  [9]="tigress "$parallel3
  [10]="tomsimilarity -m bicor")

# usage
USAGETXT=\
"
$0 <expression-matrix.tsv> <genes.tsv>

Methods to be run: ${inference[@]}

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

if [ -z "$mail"]; then
  abort "Edit the script and add your email"
fi

if [ $(wc -l $2 | cut -d" " -f1) -ne 1 ]; then
  abort "The second file should have only one line"
fi

# MIGHT ALSO consider THIS instead of the above. We could have 1 gene as a subset
# GENES=...
# HEADLESS=...
# if [[ $(stat -c "%s" $GENES) -gt $(stat -c "%s" $HEADLESS) ]]
# then
#   echo "The file size of column headers can't be greater than that of "
#   echo "expression values"
#   exit 1
#  fi

# create dirs
if [ ! -d $resultDir ]; then
  mkdir $resultDir
fi

if [ ! -d $resultDir/mi ]; then
  mkdir $resultDir/mi
fi

# Find the number of genes to set the batch size
# Rules: set -B to
# 1) the number of genes if you are using a single node
# 2) if using more nodes split the number of genes by the number of nodes
# 3) for batch mode, set it to batch
ngenes=$(cat $2 | wc -w)
default="-B $ngenes"

# the number of nodes depends on the arguments list and default

optionB=(
  [0]=$default
  [1]=$default
#  [2]="-B $(expr $(expr $ngenes '/' $(expr $nodes "*" $genie3Multiplier )) '+' 1)"
  [2]=$default
  [3]=$default
  [4]=$default
  [5]=""
  [6]=""
  [7]=$default
  [8]=""
#  [9]="-B $(expr $(expr $ngenes '/' $(expr $nodes "*" $tigressMultiplier)) '+' 1)"
  [9]=$default
  [10]="")

# Set the number of OMP threads
# default to 1
# set to -n (number of cores for pearson, spearman, and pcor
default=
ompThread=(
  [0]=$default
  [1]=$default
  [2]=$default
  [3]=$default
  [4]=$default
  [5]=$correlationNCPUs
  [6]=$correlationNCPUs
  [7]=$default
  [8]=$correlationNCPUs
  [9]=$default
  [10]=$correlationNCPUs)

# Set dependencies
# Both ARACNE and CLR calculate the RAW MI
# We will have ARACNE export it and CLR use it.
deps=(
  [1]=0
)

# Define a JobIDs array
jobIDs=()

# Create template script and submit
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}
  if [ ${run[$i]} -eq 1 ]; then
    if [ ! -f $resultDir/$inf/$inf.tsv ]; then
      mkdir -p $resultDir/$inf
      echo "#!/bin/bash" > $resultDir/$inf/$inf.sh
      if [ -f $resultDir/$inf/$inf.json ]; then
        svr="--resume-from $resultDir/$inf/$inf.json"
      else
        svr="--save-resume $resultDir/$inf/$inf.json"
      fi
      if [ -z ${ompThread[$i]} ]; then
 	      echo "unset OMP_NUM_THREADS" >> $resultDir/$inf/$inf.sh
      else
    	  echo "export OMP_NUM_THREADS=${ompThread[$i]}" >> $resultDir/$inf/$inf.sh
    	  svr=
      fi
	    echo "srun ${command[$i]} ${optionB[$i]} -i $1 -g $2 $svr -o $resultDir/$inf/$inf.tsv" >> $resultDir/$inf/$inf.sh
      
      # Handle dependencies
      dep=
      if [ ! -z ${deps[$i]} ]; then
        dep="-d afterok:${jobIDs[${deps[$i]}]}"
      fi

      dep=$(sbatch --mail-type=END,FAIL --mail-user=$mail -A $account -J $inf $dep \
      -e $resultDir/$inf/$inf.err -o $resultDir/$inf/$inf.out -p $partition ${arguments[$i]} $resultDir/$inf/$inf.sh)

      jobIDs[$i]=${dep//[^0-9]/}
    fi
  fi
done

