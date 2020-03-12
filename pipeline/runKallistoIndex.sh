#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=8GB

set -eux

module load bioinfo-tools kallisto

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <fasta file> <out dir>
    Note: The database filename defaults to the input basename
" 
    exit 1
}

## VARS
KMER=""


## get the options
while getopts k: option
do
    case "$option" in
	k) KMER="-k $OPTARG";;
	\?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

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

# get the extension
ext="${1##*.}"

if [ "$ext" == "gz" ]; then
  inxName=$2/$(basename ${1/.f*a.gz/}).inx
else
  inxName=$2/$(basename ${1/.f*a/}).inx
fi

# running
kallisto index $KMER -i $inxName $1
