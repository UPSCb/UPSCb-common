#!/bin/bash
#SBATCH -A facility
#SBATCH -t 4:00:00
#SBATCH -p core -n 2
#SBATCH --mail-type=ALL

set -ex

# variables
# we use 2 because of hyperthreading
# to force the use of a virtual core instead of a logical one, we could
# sbatch -n 1 -c 1 OMP_NUM_THREADS=1 seidr backbone -O 1
CPU=2

# helper functions
source ../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Usage: $0 <seidr file> <threshold> <output filename>
  
  Note: the threshold is the quantile value from a normal distribution,
  so a backbone of 1% is qnorm(0.99) = 2.33. 10% is 1.28. etc.
"

# sanity
isExec seidr
if [ $? -ne 0 ]; then
  abort "seidr is not available. Install it, or load the module"
fi

if [ $# -ne 3 ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

if [ ! -d $(dirname $3) ]; then
  abort "The third argument directory needs to exist"
fi

# run
export OMP_NUM_THREADS=$CPU
seidr backbone -F $2 -o $3 $1

