#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

set -ex

module load bioinfo-tools blast/2.2.26

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <fasta file> <out dir>
    Options:
            -p the type of file T/F (default to F: nucleotide; use T for protein)
            -t the db title
    Note: The database filename defaults to the input basename
" 
    exit 1
}

## VARS
OPTIONS=""
TYPE="F"

## get the options
while getopts p:t: option
do
    case "$option" in
	t) OPTIONS="$OPTIONS -t $OPTARG";;
	p) TYPE="$OPTARG";; 
	\?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

## extend the OPTIONS
OPTIONS="$OPTIONS -p $TYPE"

## we get one file and one dir as input 
if [ $# != 2 ]; then
    echo "This function takes one fasta file and one output dir as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing fasta file"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory"
    usage
fi


# running
formatdb -i $1 -o T -n $2/`basename $1` $OPTIONS
