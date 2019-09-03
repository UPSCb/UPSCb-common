#!/bin/bash -l
#SBATCH -p core
#SBATCH --mail-type=ALL

set -ex

usage(){
  echo >&2 \
"
Usage: $(basename $0) <script> <additional arguments>

This runs a script using the provided arguments. The script has to be executable.

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

# run the job
$script $@
