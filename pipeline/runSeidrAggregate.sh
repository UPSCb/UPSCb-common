#!/bin/bash
#SBATCH -A facility
#SBATCH -t 12:00:00
#SBATCH -p core -n 32
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL

set -ex

# usage
USAGETXT=\
"
  Usage: $0 <out dir> <sf file> [sf file] ... [sf file]

"
CPU=32

# sanity
if [ -z $UPSCb ]; then
  echo "Set your UPSCb environment variable"
  exit 1
fi

source $UPSCb/src/bash/functions.sh

isExec seidr

if [ $# -lt 2 ]; then
  abort "This script expects at least 2 arguments"
fi

if [ ! -d $1 ]; then
  abort "The first argument needs to be an existing directory"
fi

# run
cd $1
shift
export OMP_NUM_THREADS=$CPU
rm -f aggregated.sf
seidr aggregate -m irp -k -O $CPU $@
