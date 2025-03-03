#!/bin/bash -l
#SBATCH -n 16 -p main
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16G


## be verbose and stop on error
set -eux

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
USAGETXT=\
"
Usage: $0 [options] <input fasta> <output fasta>

Options:
   -c sequence identity threshold; default to 0.99
   -T number of threads; default to 16
   -M the memory to allocate; default to 16000 MB
   -g turn off accurate mode; default on
   -z input is compressed and output need to be compressed
Note:
   You need to set the UPSCb env. variable to your UPSCb git checkout directory
"

## defaults
IDENT=0.99
THREADS=16
MEM=16000
ACC=1
GZ=0

## options
while getopts c:T:M:gz option
do
    case "$option" in
	c) IDENT=$OPTARG;; 
	T) THREADS=$OPTARG;;
	M) MEM=$OPTARG;;
	g) ACC=0;;
  z) GZ=1;;
	\?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## arguments
if [ $# != 2 ]; then
    abort "This function takes 2 arguments: the input fasta file and the output filename"
fi

if [ ! -f $1 ]; then
    abort "The first argument must be a valid fasta file"
fi

if [ ! -d $(dirname $2)]; then
    abort "The second argument must be an output filename which parent directory must exist"
fi

## tool
isExec cd-hit-est

if [ "$GZ" == "1" ]; then
  ## create temporary file
  tmp=$(tempfile)

  ## decompress
  gunzip -c $1 > $tmp
else
  tmp=$1
fi

## command
# -d 0 is to ensure that the whole ID is kept in the cluster file
cd-hit-est -i $tmp -o $2 -c $IDENT -T $THREADS -M $MEM -g $ACC -d 0 -p 1

if [ "$GZ" == "1" ]; then
  gzip $2
  rm $tmp
fi
