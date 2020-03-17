#!/bin/bash
#SBATCH -A facility
#SBATCH -t 4:00:00
#SBATCH -p core -n 1
#SBATCH --mem=16GB

set -ex

# Options
ALL="-a"
INDEX=
POINTS="-p 1000"

# usage
USAGETXT=\
"
  Usage: $0 [options] <seidr file> <positive-gold-standard> <negative-gold-standard> <output filename>
  
  Options:
        -a (default, set the flag to unset) process all algorithms
        -i index of the algorith to proceed
        -p the number of points to keep (default: $POINTS)
"

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

isExec seidr

# Get the options
while getopts ai:p: option
do
    case "$option" in
        a) ALL="";;
        i) INDEX="-i $OPTARG";;
        p) POINTS="-p $OPTARG";;
        \?) ## unknown flag
		    abort;;
    esac
done
shift `expr $OPTIND - 1`

OPTIONS="$ALL $INDEX $POINTS"

if [ $# -ne 4 ]; then
  abort "This script expects 4 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing file"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be an existing file"
fi

if [ ! -f $3 ]; then
  abort "The third argument needs to be an existing file"
fi

if [ ! -d $(dirname $4) ]; then
  abort "The fourth argument directory needs to exist"
fi

# run
seidr roc -n $1 -g $2 -x $3 $OPTIONS > $4
