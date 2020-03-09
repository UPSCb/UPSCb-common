#!/bin/bash

# safeguards
set -ex

# project vars
account=facility
mail=iryna.shutava@umu.se

# source helpers
source functions.sh

# Variables
inference=(aracne clr elnet genie3 llr-ensemble mi narromi pcor pearson plsnet spearman tigress)

# additional parameters (elnet is done iteratively, so the format is not the expected one: a matrix, rather an edge list
arguments=([2]="-f el" [4]="-o results/llr-ensemble/llr-ensemble.sf")
CPUs=4
MEM=64G
Time=1-00:00:00
ShortTime=1:00:00

# usage
USAGETXT=\
"
$0 <genes.tsv>

Methods to be imported from: ${inference[@]}

"

# input: gene and data
if [ $# -ne 1 ]; then
  abort "This script expects 1 argument"
fi

if [ ! -f $1 ]; then
  abort "The argument needs to be the gene names tab delimited file"
fi

# Create template script and submit
jobID=
len=${#inference[@]}
for ((i=0;i<len;i++)); do
  inf=${inference[$i]}
  if [ "$inf" == "elnet" ] && [ -d results/el-ensemble ]; then
    # create a job to cat el-ensemble into elnet
    if [ ! -d results/$inf ]; then
	    mkdir -p results/$inf
      echo "#!/bin/bash" > results/$inf/cat.sh
      echo "cat results/el-ensemble/*.tsv > results/$inf/$inf.tsv ">> results/$inf/cat.sh
      dep=$(sbatch -t $ShortTime -c 1  --mail-type=ALL --mail-user=$mail -A $account -J $inf-cat \
      -e results/$inf/$inf-$i.err -o results/$inf/$inf-$i.out \
      results/$inf/cat.sh)
      jobID="-d afterok:${dep//[^0-9]/}"
    fi
  fi
  if [ ! -f results/$inf/$inf.sf ]; then
    if [ -f results/$inf/$inf.tsv ]; then

      ./generate_import_script.py \
      -i results/$inf/$inf.tsv -g $1 -c $CPUs ${arguments[$i]} > results/$inf/${inf}-import.sh

      sbatch -t $Time --mail-type=ALL --mail-user=$mail -A $account -J import-$inf $jobID \
	    -e results/$inf/${inf}-import.err -o results/$inf/${inf}-import.out \
	    -c $CPUs --mem=$MEM results/$inf/${inf}-import.sh
    else
      if [ "$jobID" != "" ]; then
        # most likely the tsv file for elnet won't exist, so give it a chance (as there is a dep)
	      # we touch the file so the generate import script does not fail
	      touch results/$inf/$inf.tsv

        ./generate_import_script.py \
        -i results/$inf/$inf.tsv -g $1 -c $CPUs ${arguments[$i]} > results/$inf/${inf}-import.sh

        sbatch -t $Time --mail-type=ALL --mail-user=$mail -A $account -J import-$inf $jobID \
	      -e results/$inf/${inf}-import.err -o results/$inf/${inf}-import.out \
		-c $CPUs --mem=$MEM results/$inf/${inf}-import.sh
      else
        echo "There is no tsv file for $inf"
      fi
    fi
    # ${arguments[$i]}
  fi
done

