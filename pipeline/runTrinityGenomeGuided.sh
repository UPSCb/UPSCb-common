#!/bin/bash -l
#SBATCH -p core
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
##SBATCH --mem=128G
## time too for large files
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL

## stop on error
set -eu

## check the options if any
PROC=20
INTRON=11000
MEM="128G"
STSPEC=
LONGREADS=
NORM=

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
USAGETXT=\
"
	Usage: runTrinityGenomeGuided.sh [options] <singularity container> <out dir> <BAM file>

	Options:
	              -i max intron size
                -l long reads (pac bio or other) as a fasta file
                -m mem requirement (in GB, default 128)
                -n no digital read normalisation
		            -p number of threads to use (default 16)
		            -r strand specific (forward first (ISF as per salmon notation))
                -s strand specific (Illumina dUTP protocol)

  Note:
    -r and -s are mutually exclusive; they overload the same parameter; this is not checked for.
"

## get the options
while getopts l:m:np:rsi: option
do
    case "$option" in
      i) INTRON=$OPTARG;; 
      l) LONGREADS="$OPTARG";;
      n) NORM="--no_normalize_reads";;
      m) MEM=$OPTARG;;
      p) PROC=$OPTARG;;
      r) STSPEC="--SS_lib_type FR";;
      s) STSPEC="--SS_lib_type RF";;
	 	\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

# sanity
## we get one dir and one file
[[ $# != 3 ]] && abort "This function needs 3 arguments"

[[ ! -f $1 ]] && abort "The first argument (singularity container) needs to be an existing file"

[[ ! -d $2 ]] && abort "The second argument (output dir) needs to be an existing directory"

[[ ! -f $2 ]] && abort "The third argument (alignment BAM) needs to be an existing file"

## check args
if [ ! -z $LONGREADS ]; then
	[[ ! -f $LONGREADS ]] && abort "-l should point to an existing fasta file"
	LONGREADS="--long_reads $LONGREADS"
fi

## run trinity
singularity exec $1 \
Trinity --genome_guided_bam $3 --genome_guided_max_intron $INTRON \
--max_memory $MEM --CPU $PROC --output $2 $LONGREADS $NORM $STSPEC

