#!/bin/bash -l
#SBATCH -p core
#SBATCH --mail-type=ALL

set -ex

usage(){
  echo >&2 \
"
Usage: $(basename $0) <script> <file list> <additional arguments>

This runs an array of the script using as arguments:
(i) iteratively a line of file list (i.e. specific arguments, such as filename(s))
(ii) common arguments provided as additional arguments on the command line


"
exit 1
}

# check that the script exists
if [ ! -f $1 ]; then
  echo 
  usage
fi
script=$1
shift

# check that the file list exists
if [ ! -f $1 ]; then
  echo 
  usage
fi

# read the file list
readarray -t array < $1
shift

# run the jobs
bash $script ${array[$SLURM_ARRAY_TASK_ID]} $@
