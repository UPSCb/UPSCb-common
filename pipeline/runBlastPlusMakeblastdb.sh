#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

set -ex

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <fasta file> <out dir>
    Options:
            -p the type of file nucl/prot (default to nucl)
            -t the db title
    Note: The database filename defaults to the input basename
" 
    exit 1
}

## VARS
OPTIONS=""
TYPE="nucl"

## get the options
while getopts p:t: option
do
    case "$option" in
	p) 
	    case $OPTARG in
		nucl)TYPE="$OPTARG";;
		prot)TYPE="$OPTARG";;
		*)
		    echo "You provided an unsupported file type, use either 'nucl' or 'prot'";
		    usage;;
	    esac;;
	t) OPTIONS="-title $OPTARG";;
	\?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

## extend the options
OPTIONS="$OPTIONS -dbtype $TYPE"

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
makeblastdb -in $1 -out $2/`basename $1` $OPTIONS
