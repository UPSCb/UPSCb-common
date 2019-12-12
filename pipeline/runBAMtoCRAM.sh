#!/bin/bash
#SBATCH -p core -n 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

set -ex

# usage txt
export USAGETXT=\
"
	Usage: $0 [options] <genome fasta> <bam file>
"

# common function
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# options
while getopts h option
do
  case "$option" in
      h) usage;;
      ?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# check the arguments
if [ "$#" != 2 ]; then
  abort "This function expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "This function expects a fasta file as first argument"
fi

if [ ! -f $2 ]; then
  abort "This function expects a bam file"
fi

# check tool
isExec samtools

# run
in=$2
fasta=$1
out=${in/.bam/.cram}

samtools view -C -T $fasta -o $out $in

samtools index $out

if [ -f $out ]; then
    rm $in
fi

if [ -f $in.bai ]; then
    rm $in.bai
fi
