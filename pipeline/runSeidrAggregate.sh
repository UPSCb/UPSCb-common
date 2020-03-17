#!/bin/bash
#SBATCH -A facility
#SBATCH -t 12:00:00
#SBATCH -p core -n 32
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL

set -ex
DIRECTIONALITY="-k"
FORCE=
METHOD="-m irp"

# usage
USAGETXT=\
"
  Usage: $0 [options] <out dir> <sf file> [sf file] ... [sf file]
  
  Options:
          -f  force overwrite output
          -k  keep the directionality, default: true
          -m  method, default: irp
          
"
CPU=32

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

isExec seidr

# Get the options
while getopts fkm: option
do
    case "$option" in
        f) FORCE="-f";;
        k) DIRECTIONALITY=;;
        k) METHOD="-m $OPTARG";;
        \?) ## unknown flag
		    abort;;
    esac
done
shift `expr $OPTIND - 1`

OPTIONS="$FORCE $DIRECTIONALITY $METHOD"

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
#rm -f aggregated.sf

seidr aggregate $OPTIONS -O $CPU $@
