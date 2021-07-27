#!/bin/bash -l
#SBATCH -A SNIC2021-5-200
#SBATCH -t 12:00:00
#SBATCH -n 28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicolas.delhomme@slu.se
#SBATCH -o aggregate.out
#SBATCH -e aggregate.err

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

CPU=28

# source
source functions.sh

# modules
EXEC=/pfs/proj/nobackup/fs/projnb10/snic2019-35-44/software/seidr/build
source $EXEC/sourcefile

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
