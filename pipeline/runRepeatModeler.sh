#!/bin/bash -l
#SBATCH -J repeatModeler
#SBATCH -p core
#SBATCH -c 8
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x
CPU=8

# usage txt
export USAGETXT=\
"
	Usage: $0 [options] <repeatModelerDB> <outdir>
	
	Options: -c the number of CPU to use
"

# common function
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# options
while getopts c:h option
do
  case "$option" in
      c) CPU=$OPTARG;;
      h) usage;;
      ?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# check the arguments
if [ "$#" != 2 ]; then
  abort "This function expects at least two arguments"
fi

genome=$1
shift
if [ ! -f $genome.nhr ]; then
  abort "This function expects a database prefix as first argument"
fi

outdir=$1
shift
if [ ! -d $outdir ]; then
  abort "The output directory needs to exist"
fi

# run
cd $outdir
RepeatModeler -pa $(expr $CPU - 1) -database $genome
