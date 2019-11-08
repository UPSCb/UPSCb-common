#!/bin/bash -l
#SBATCH -J repeatmasker
#SBATCH -p core
#SBATCH -c 8
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# usage txt
export USAGETXT=\
"
	Usage: $0 <genome> <engine> <outdir> [options to RM]
  Note: a default option is -qq, set another option to overwrite
"

OPT="-qq"

# common function
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

if [ "$#" -lt 3 ]; then
  abort "This function expects at least three arguments"
fi

genome=$1
shift
if [ ! -f $genome ]; then
  abort "This function expects a fasta file as first argument"
fi

engine=$1
shift

outdir=$1
shift
if [ ! -d $outdir ]; then
  abort "The output directory needs to exist"
fi

if [ "$#" -gt 1 ]; then
  OPT=$@
fi

# run
RepeatMasker $genome -e $engine -pa 8 -dir $outdir $OPT


