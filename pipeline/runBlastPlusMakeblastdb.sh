#!/bin/bash -l
#SBATCH -p main -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

set -eux

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# a usage function
USAGETXT=\
"
    Usage: $0 [options] <singularity blast container> <fasta file> <out dir>
    Options:
            -p the type of file nucl/prot (default to nucl)
            -t the db title
    Note: The database filename defaults to the input basename
"

# VARS
OPTIONS=""
TYPE="nucl"

# get the options
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
	\?) # unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

# extend the options
OPTIONS="$OPTIONS -dbtype $TYPE"

# We rely on singularity

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# we get one file and one dir as input 
[[ $# != 3 ]] && abort "This function takes a singularity container, one fasta file and one output dir as argument"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity container"

[[ ! -f $2 ]] && abort "The second argument needs to be an existing fasta file"

[[ ! -d $3 ]] && abort "The third argument needs to be an existing directory"
    
# Check if the input is compressed (based on file extension)
if [ "${2##*.}" == "gz" ]; then
 in=$(mktemp)
 gunzip -c $2 > $in
else
 in=$2
fi

# running
singularity exec $1 makeblastdb -in $in -out $3/`basename $2` $OPTIONS
