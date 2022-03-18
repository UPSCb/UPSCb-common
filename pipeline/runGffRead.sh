#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t -00:30:00
#SBATCH --mail-type=ALL

# stop on error, be verbose
set -eu

## defaults
TYPE="gtf"

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# a usage function
USAGETXT= \
"
    Usage: $0 [options] <gff/gtf file>

    Options:
            -t the output type: gtf or gff [default to $TYPE]

    Note: 
         1) The tool is meant to convert gtf to gff and vice-versa. Hence to output 'gtf' a 'gff/gff3' input is expected.
         2) The output filename is the input filename with the converted extension. If the file exists, it is not overwritten
            but an additional '.new' extension is added.
" 

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
[[ $# != 2 ]] && abort "This function takes two arguments, the singularity gffread container and the gff/gtf file path"

## check that the files exist
[[ ! -f $1 ]] && abort "The argument shoud be a valid file path to a singularity container"
[[ ! -f $2 ]] && abort "The argument shoud be a valid file path to a gtf/gff file"

## enforce singularity
[[ -z $SINGULARITY_BINDPATH ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

## check the option
## and define the output
extension="${2##*.}"
option=
case $TYPE in
    gtf) 
	    [[ "$extension" != "gff" ]] && [[ "$extension" != "gff3" ]] && \
	    abort "The expected extension for your input file is 'gff' or 'gff3'"
    	option="-T";;    
    gff) 
	    [[ "$extension" != "gtf" ]] && abort "The expected extension for your input file is 'gtf'";;
    *) abort "The type must be one of 'gtf' or 'gff'"
esac

# define the output and check for existence
# we do not overwrite
out=${2//.$extension/.$TYPE}
[[ -f $out ]] && echo "The output file $out already exists, renaming the new file to $out.new" && out=$out.new

## run it
singularity exec $1 gffread $2 -o $out $option $2
