#!/bin/bash
#SBATCH -p node -n 20
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

set -ex

# usage txt
export USAGETXT=\
"
	Usage: $0 [options] <target fasta> <query fasta> <outdir>
	
	Options: -c the number of CPU to use
"

# common function
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# vars
CPU=20

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
if [ "$#" != 3 ]; then
  abort "This function expects three arguments"
fi

if [ ! -f $1 ]; then
  abort "This function expects a fasta file as first argument"
fi

if [ ! -f $2 ]; then
  abort "This function expects a fasta file as second argument"
fi

if [ ! -d $3 ]; then
  abort "The output directory needs to exist"
fi

isExec minimap2

# run
cd $outdir
minimap2 -t $CPU -x map-pb -a $1 $2 > $3/$(basename ${1/.fa*/})
