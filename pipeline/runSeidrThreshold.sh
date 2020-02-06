#!/bin/bash
#SBATCH -A facility
#SBATCH -p core -n 16
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=96GB

# variables
CPU=16

# usage
USAGETXT=\
"
  Usage: $0 <seidr file> <output filename>
  
"

# sanity
if [ -z $UPSCb ]; then
  echo "Set your UPSCb environment variable"
  exit 1
fi

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

isExec seidr

if [ $# -ne 2 ]; then
  abort "This script expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

if [ ! -d $(dirname $2) ]; then
  abort "The second argument directory needs to exist"
fi

# run
export OMP_NUM_THREADS=$CPU
seidr threshold -n 10000 -m 0 -M 1 -O $CPU -o $2 $1
