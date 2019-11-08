#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 12
## no mail at the time
###SBATCH --mail-type=ALL

## safeguards
set -eux

## functions
source ../UPSCb-common/src/bash/functions.sh

## options default
FMT=16
EVALUE=1e-5
PROC=12
OPTIONS=""

## usage
USAGETXT=\
"
	Usage: $0 [options] <blast command> <fasta file> <index> <out dir>

	Options:
    -e the e-value (default to $EVALUE)
		-p number of threads to use (default $PROC)
		-b blast options (example -task blastn instead of megablast)

  Note:
	     The default format is 16 (single file in XML2)
"

## get the options
while getopts e:p:b: option
do
    case "$option" in
	    e) EVALUE=$OPTARG;;
	    p) PROC=$OPTARG;;
	    b) OPTIONS=$OPTARG;;
		  \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    abort "This function takes one blast command, one file, one index and one directory as arguments"
fi

BLAST=$1

# Executable
isExec $BLAST
if [ $? -ne 0 ]; then
  abort "$1 is not available. Install it, or load the module"
fi

# Extension
EXT="nhr"
case "$1" in
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
	abort "Unknown blast command. Not one of blastn, blastp, blastx, tblastn, tblastx. Aborting."
esac

if [ ! -f $2 ]; then
    abort "The second argument needs to be the fasta file"
fi

if [ ! -f $3.$EXT ] && [ ! -f $3.${EXT//hr/al} ]; then
    abort "The third argument needs to be the basename of the BLAST index - including the path"
fi

if [ ! -d $4 ]; then
    abort "The forth argument needs to be the output directory"
fi

## run BLAST
echo Blasting
fnam=$(basename $3)_$(basename ${2//.f*a*/.xml})
$1 -db $3 -query $2 -out $4/$fnam -evalue $EVALUE $OPTIONS -num_threads $PROC -outfmt $FMT

##
echo Done


