#!/bin/bash -l
#SBATCH -p all
#SBATCH -n 1
## no mail at the time
###SBATCH --mail-type=ALL

## stop on error
set -e

## be verbose and extend the commands
set -x

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## check the options if any
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
EVALUE=1e-5
PROC=1
OPTIONS=""

## usage
USAGETXT=\
"
	Usage: $0 [options] <singularity container><blast command> <fasta file> <index> <out dir>

	Options:
                -f the output format (default to $FMT)
                -e the e-value (default to $EVALUE)
		-p number of threads to use (default $PROC)
		-b blast options (example -task blastn instead of megablast)

        Note:
	     Providing the -f option is broken at the moment
"

## get the options
while getopts e:f:p:b: option
do
        case "$option" in
	    f) FMT=$OPTARG;;
	    e) EVALUE=$OPTARG;;
	    p) PROC=$OPTARG;;
	    b) OPTIONS=$OPTARG;;
		\?) ## unknown flag
		abort "unknown option";;
        esac
done
shift `expr $OPTIND - 1`

## we get one container, two dir and two files as input
[[ $# != 5 ]] && abort "This function takes one container, one blast command, one file, one index and one directory as arguments"

[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing ncbi-blast container"

BLAST=
EXT="nhr"
case "$2" in
    blastn)
	;;
    blastp)
	EXT="phr"
	;;
    blastx)
	EXT="phr"
	;;
    tblastn)
	;;
    tblastx)
	;;
    *)
	echo "Unknown blast command. Not one of blastn, blastp, blastx, tblastn, tblastx. Aborting."
	usage
esac
BLAST=$2

[[ ! -f $3 ]] && abort "The third argument needs to be the fasta file"

[[ ! -f $4.$EXT ]] && [[ ! -f $4.${EXT//hr/al} ]] && abort "The fourth argument needs to be the basename of the BLAST index - including the path"

[[ ! -d $5 ]] && abort "The fifth argument needs to be the output directory"

## run BLAST
echo Blasting
fnam=$(basename $4)_$(basename ${3//.f*a*/.blt})
singularity exec $1 $2 -db $4 -query $3 -out $5/$fnam -evalue $EVALUE $OPTIONS -num_threads $PROC -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

##
echo Done





