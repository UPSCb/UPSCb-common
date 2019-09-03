#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL

## stop on error
set -e

## modules
module load bioinfo-tools Picard-tools

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <fasta file>

    Options:
            -s the species to be added to the SP tag
" 
    exit 1
}

## VARS
OPTIONS=""

## get the options
while getopts s: option
do
    case "$option" in
	s) OPTIONS="SPECIES=$OPTARG";;
	\?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`


## we get one file as input
if [ $# != 1 ]; then
    echo "This function takes one fasta file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing fasta file"
    usage
fi

## create the output
out=${1//.f*a/.dict}

## create the index
java -jar  $PICARD_TOOLS_DIR/picard.jar CreateSequenceDictionary REFERENCE=$1 OUTPUT=$out $OPTIONS
