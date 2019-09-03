#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## usage
usage(){
echo >&2 \
"
	Usage: runHTSeq.sh [options] <out dir> <in.bam> <in.gff>	

	Options:
                -i precise the IDATTR
                   default to 'Parent', but e.g. should be 'pacid' 
                   for the P. trichocarpa gene exon gff3 file
                -s is the protocol stranded?
                   default to FALSE
                -a are we counting antisense transcripts?
                    default to FALSE, only active in combination with -s
                -t Chose attribute to count in the gff3 file default is exon
        Note:
                BAM file are expected to be sorted by position
                Only HTSeq 0.6+ version(s) are supported
"
	exit 1
}

echo Loading modules
module load bioinfo-tools htseq

## check the version
#if [ `htseq-count --help | grep -c "version 0.6"` -ne 1 ]; then
#    echo Only HTSeq version 0.6+ are supported
#    usage
#fi

## options
IDATTR="Parent"
stranded=0
antisense=0
t="exon"

## get the options
while getopts ai:st: option
do
        case "$option" in
      a) antisense=1;;
	    i) IDATTR=$OPTARG;;
	    s) stranded=1;;
      t) t=$OPTARG;;
	    \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 3 ]; then
    echo "This function takes one directory, one bam and one gff3 file as arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing bam file"
    usage
fi
nam=`basename ${2//.bam/}`

if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing gff3 file"
    usage
fi

if [ $t == "CDS" ]; then
  echo "Warning: the CDS option require the CDS feature to be capital in you gff3 file"
fi

## get the count table
if [ $stranded == 0 ]; then
  if [ $antisense == 1 ]; then
    echo "The antisense only works in conjunction with the -s option" >&2
  fi
## since we are not using strand specific, go for the union
    htseq-count -f bam -r pos -m union -s no -t $t -i $IDATTR $2 $3 > $1/$nam.txt
else
  ## normal counting
  if [ $antisense == 0 ]; then
    htseq-count -f bam -r pos -m intersection-nonempty -s reverse -t $t -i $IDATTR $2 $3 > $1/$nam.txt
  else
    htseq-count -f bam -r pos -m intersection-nonempty -s yes -t $t -i $IDATTR $2 $3 > $1/$nam.txt
  fi
fi

