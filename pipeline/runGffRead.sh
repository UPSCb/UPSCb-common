#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t -00:30:00
#SBATCH --mail-type=ALL

## stop on error, be verbose
set -ex

## load modules
module load bioinfo-tools cufflinks

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <gff/gtf file>

    Options:
            -t the output type: gtf or gff [default to gtf]

    Note: 
         1) The tool is meant to convert gtf to gff and vice-versa. Hence to output 'gtf' a 'gff/gff3' input is expected.
         2) The output filename is the input filename with the converted extension. If the file exists, it is not overwritten
            but an additional '.new' extension is added.
" 
    exit 1
}

## defaults
TYPE="gtf"

## get the options
while getopts t: option
do
    case "$option" in
	t) TYPE=$OPTARG;;
	\?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

## check the number of arguments
if [ $# != 1 ]; then
    echo "This function take one argument, the gff/gtf file path"
    usage
fi

## check that the file exists
if [ ! -f $1 ]; then
    echo "The argument shoud be a valid file path to a gtf/gff file"
    usage
fi

## check the option
## and define the output
extension="${1##*.}"
option=
case $TYPE in
    gtf) 
	if [ "$extension" != "gff" ] && [ "$extension" != "gff3" ]; then
	    echo "The expected extension for your input file is 'gff' or 'gff3'"
	    usage
	fi
	option="-T"
	;;    
    gff) 
	if [ "$extension" != "gtf" ]; then
	    echo "The expected extension for your input file is 'gtf'"
	    usage
	fi
	;;
    *) echo "The type must be one of 'gtf' or 'gff'"
	usage;;
esac

# define the output and check for existence
# we do not overwrite
out=${1//.$extension/.$TYPE}
if [ -f $out ]; then
    echo "The output file $out already exists, renaming the new file to $out.new"
    out=$out.new
fi 

## run it
gffread $1 -o $out $option $1
