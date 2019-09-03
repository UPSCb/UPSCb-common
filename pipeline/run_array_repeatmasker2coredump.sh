#!/bin/bash -l

#SBATCH -J arrayrepeatmasker.job
#SBATCH -p core
#SBATCH -c 8
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools 
module load RepeatMasker

# usage function

usage(){
echo >&2 \
"
	Usage: $0 <link to genome>

"
	exit 1
}

# run the command # -dir $3 
#RepeatMasker $1.$SLURM_ARRAY_TASK_ID.masked_core.fasta -e $2 -pa $3 -qq -lib $4
RepeatMasker $1 -e $2 -pa $3 -qq -lib $4



