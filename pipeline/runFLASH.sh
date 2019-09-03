#!/bin/bash -l
#SBATCH -p node 
#SBATCH -n 16
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=ALL

## abort on error
set -e

## usage
usage(){
echo >&2 \
"Usage: 
    $0 [option] <fwd fastq file> <rev fastq file> <output dir> [additional arguments]

Options:
    -r      read length (default 151)
    -f      fragment size (default 200)
    -s      fragment size standard deviation (default 25)
    -t      number of threads
    -c      merge uncombined read 1 with combined
    -o      output fasta

Notes:
    The -o option requires the fastx toolkit to be available
"
    exit 1
}

## options
frag=250
fragsd=25
read=151
thread=16
combine=0
fasta=0

while getopts "cof:r:s:t:" opt; do
  case $opt in
    c) combine=1;;
    o) fasta=1;;
    f) frag=$OPTARG;;
	  r) read=$OPTARG;;
	  s) fragsd=$OPTARG;;
    t) thread=$OPTARG;;
    \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

## check file 1
if [ ! -f $1 ]; then
    echo "The first argument must be the valid file name of the forward fastq file."
    usage
fi
f=$1
shift

## check file 2
if [ ! -f $1 ]; then
    echo "The second argument must be the valid file name of the reverse fastq file."
    usage
fi
r=$1
shift

## check dir
if [ ! -d $1 ]; then
    echo "The third argument must be an existing directory."
    usage
fi
d=$1
shift

## get the prefix
fnam=`basename ${f//_[1,2].f*q.gz/}`

## output
flash -r $read -f $frag -s $fragsd -z -t $thread -o $fnam -d $d $f $r $@

if [ $combine -eq 1 ]; then
  cat $d/${fnam}.notCombined_1.fastq.gz $d/${fnam}.extendedFrags.fastq.gz > $d/${fnam}_flash.fq.gz
fi

if [ $fasta -eq 1 ]; then
  if [ $combine -eq 1 ]; then
    gunzip -c $d/${fnam}_flash.fq.gz | fastq_to_fasta -o $d/${fnam}_flash.fa
  else
    gunzip -c  $d/${fnam}.extendedFrags.fastq.gz | fastq_to_fasta -o  $d/${fnam}.extendedFrags.fa
  fi
fi
