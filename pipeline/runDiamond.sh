#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 20
## no mail at the time
#SBATCH --mail-type=ALL

## stop on error, unset and be verbose
set -eux

## load the modules
module load bioinfo-tools diamond

# check for the env var
UPSCb=${UPSCb:-}
if [ -z $UPSCb ]; then
    echo "You need to set your UPSCb environment variable to your UPSCb checkout directory"
    exit 1
fi

## source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## check the options if any
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids"
PROC=20
OPTIONS=""

## usage
USAGETXT=\
"
	Usage: $0 [options] <diamond command> <fasta file> <index> <out dir>

	Options:
	        -i taxid file (from NCBI taxonomy: prot.accession2taxid.gz)
		    -m more sensitive
		    -n nodes.dmp file (from NCBI taxonomy)
		    -o additional diamond options
		    -p number of threads to use (default $PROC)
		    -s include all titles
		    -t names.dmp  file (from NCBI taxonomy)
		    
    Note:
            the output format defaults to $FMT
"

## get the options
while getopts i:mn:o:p:st: option
do
        case "$option" in
        i) OPTIONS="$OPTIONS --taxonmap $OPTARG";;
        m) OPTIONS="$OPTIONS --more-sensitive";;
        n) OPTIONS="$OPTIONS --taxonnodes $OPTARG";;
	    p) PROC=$OPTARG;;
	    o) OPTIONS="$OPTIONS $OPTARG";;
	    s) OPTIONS="$OPTIONS --salltitles";;
	    t) OPTIONS="$OPTIONS --taxonnames $OPTARG";;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function takes one blast command, one file, one index and one directory as arguments"
    usage
fi

CMD=(blastp blastx)
if [ $(containsElement $1 "${CMD[@]}") -eq 1 ]; then
  abort "Unknown command: $1, it should be one of ${CMD[@]}"
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be the fasta file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be the diamond index"
    usage
fi

if [ ! -d $4 ]; then
    echo "The forth argument needs to be the output directory"
    usage
fi

## run BLAST
echo "Blasting but faster"
fnam=$(basename $3)_$(basename ${2//.f*a*/.blt})

diamond $1 -d $3 -p $PROC -q $2 -o $4/$fnam -f $FMT $OPTIONS

##
echo Done


