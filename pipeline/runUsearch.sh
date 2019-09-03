#!/bin/bash -l

## be verbose and stop on error
set -ex

## usage
usage(){
echo >&2 \
"
Usage: $0 [options] <input fasta> <output prefix>

Options:
   -c sequence identity threshold; default to 0.99
   -T number of threads; default to 16
   -i the idprefix size; default to 0
Note:
   You need to set the UPSCb env. variable to your UPSCb git checkout directory
"
exit 1
}

## defaults
IDENT=0.99
THREADS=16
IDP=0

## options
while getopts c:T:i: option
do
    case "$option" in
	c) IDENT=$OPTARG;; 
	T) THREADS=$OPTARG;;
	i) IDP=$OPTARG;;
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

## command
usearch -cluster_fast $1 -id $IDENT -uc $2.uc -idprefix $IDP --centroids $2.fasta -threads $THREADS
