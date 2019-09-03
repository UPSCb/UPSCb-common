#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL

## set error and out
set -e -x

## module load
module load python/2.7.1

usage(){
echo >&2 \
"
	Usage: runBESST-RNA.sh [option] <genome fa file> <bam file> <out dir>
	
	Options:
                -e minimum support per link (3)
                -k minimum scaffold size (200)
		-m minimum aln quality (10)
		-t intron size (20000)
                -z max coverage (30)
	
	Notes:
		The BESST_RNA env. var. needs to point to your BESST_RNA Git checkout dir.
"
	exit 1
}

## the defaults
EDGE=3
MINSIZE=200
MINQUAL=10
INTRONSIZE=20000
MAXCOV=30

## get the options
while getopts e:k:m:t:z: option
do
        case "$option" in
	    e) EDGE=$OPTARG;;
	    k) MINSIZE=$OPTARG;;
            m) MINQUAL=$OPTARG;;
            t) INTRONSIZE=$OPTARG;;
	    z) MAXCOV=$OPTARG;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## check the arguments
if [ $# != 3 ] ; then
    echo "This script takes 3 arguments"
    usage
fi

if [ ! -f $1 ] ; then
    echo "The first argument should be the genome fasta file."
    usage
fi

if [ ! -f $2 ] ; then
    echo "The second argument should be the BAM alignment file."
    usage
fi

if [ ! -d $3 ] ; then
    echo "The third argument should be the output directory"
    usage
fi

if [ -z $BESST_RNA] ; then
    echo "THE BESST_RNA env. var. is not set"
    usage
fi

## create a BESST-RNA runner with a mapq option
python $BESST_RNA/src/Main.py 1 -c $1 -f $2 -o $3 -e $EDGE -T $INTRONSIZE -k $MINSIZE -d 1 -g 0 -z $MAXCOV --mapq $MINQUAL
