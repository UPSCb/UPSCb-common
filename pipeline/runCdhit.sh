#!/bin/bash -l
#SBATCH -n 16 -p core
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16G


## be verbose and stop on error
set -ex

## usage
usage(){
echo >&2 \
"
Usage: $0 [options] <input fasta> <output fasta>

Options:
   -c sequence identity threshold; default to 0.99
   -T number of threads; default to 16
   -M the memory to allocate; default to 16000 MB
   -g turn off accurate mode; default on
   -t tolerance for redundance; default 0
   -z input is compressed and output need to be compressed
Note:
   You need to set the UPSCb env. variable to your UPSCb git checkout directory
"
exit 1
}

## defaults
IDENT=0.99
THREADS=16
MEM=16000
ACC=1
TOL=0
GZ=0

## options
while getopts c:T:M:gt:z option
do
    case "$option" in
	c) IDENT=$OPTARG;; 
	T) THREADS=$OPTARG;;
	M) MEM=$OPTARG;;
	g) ACC=0;;
	t) TOL=$OPTARG;;
  z) GZ=1;;
	\?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## arguments
if [ $# != 2 ]; then
    echo "This function takes 2 arguments: the input fasta file and the output filename"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument must be a valid fasta file"
    usage
fi

if [ ! -d `dirname $2` ]; then
    echo "The second argument must be an output filename which parent directory must exist"
    usage
fi

if [ "$GZ" == "1" ]; then
  ## create temporary file
  tmp=`tempfile`

  ## decompress
  gunzip -c $1 > $tmp
else
  tmp=$1
fi

## command
# -d 0 is to ensure that the whole ID is kept in the cluster file
cd-hit -i $tmp -o $2 -c $IDENT -T $THREADS -M $MEM -g $ACC -t $TOL -d 0 -p 1

if [ "$GZ" == "1" ]; then
  gzip $2
  rm $tmp
fi
