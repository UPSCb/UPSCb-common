#!/bin/bash -l
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
## time too for large files
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

## stop on error
set -eu

## be verbose and extend the commands
set -x

## check the options if any
MEM=40G
KMER=200
PROC=20
SINGLE=0

## source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
USAGETXT=\
"
	Usage: runDigitalNormalization.sh [options] <singularity container> <left fq> <right fq> <out dir> [trinity options]
	runDigitalNormalization.sh [options] -s <singularity container> <single fq> <out dir> [trinity options]

	Options:
                -k max kmer cov (default 200)
                -m mem requirement (default 40G) (rule of thumb 1.5G/1M read pairs)
		            -p number of threads to use (default 16)
                -s for single (SE) files

        Note:
             <left fq> and <right fq> could also be files containing a list of filenames (one per line)
"

## get the options
while getopts k:m:p:s option
do
        case "$option" in
	    k) KMER=$OPTARG;;
	    m) MEM=$OPTARG;;
	    p) PROC=$OPTARG;;
	    s) SINGLE=1;;
		\?) ## unknown flag
		abort "Unknown option";;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $SINGLE -eq 0 ]; then
  [[ $# -lt 4 ]] && abort "This function takes three files and one directory as arguments for PE data"
else 
  [[ $# -lt 3 ]] && abort "This function takes two files and one directory as arguments for SE data"
fi

singularity=$1
shift
[[ ! -f $singularity ]] && abort "The first argument needs to be the trinity singularity container file"

## enforce singularity
[[ -z $SINGULARITY_BINDPATH ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

fwd=$1
shift
[[ ! -f $fwd ]] && abort "The second argument needs to be the forward (left) fastq file"

if [ $SINGLE -eq 0 ]; then
    rev=$1
    shift
    [[ ! -f $rev ]] && abort "The third argument needs to be the reverse (right) fastq file"
else
    dir=$1
    shift
    [[ ! -d $dir ]] && abort "The third argument (output dir) needs to be an existing directory"
fi

if [ $SINGLE -eq 0 ]; then
    dir=$1
    shift
    [[ ! -d $dir ]] &&abort  "The forth argument (output dir) needs to be an existing directory"
fi

## do we have file lists
LIST=""
if [ $(file --mime-type "$fwd" | grep -c "gzip$") -eq 1 ]; then
    echo "gzip formatted input"
else
	firstLine=$(head -1 $fwd)

	# TODO check for file separators :
	# if not check all files
	# if yes split and check all files
	# list and single, just concatenate
	
	if [ -f $firstLine ]; then
	  LIST="_list"
	  if [ $SINGLE -ne 0 ]; then
		  # we concatenate instead
		  for f in $(cat $fwd); do
			  l=${l:+$l,}$f
		  done
		  fwd=$l
	  fi
	fi
fi

## run trinity
if [ $SINGLE -eq 0 ]; then
	singularity exec $singularity \
	/usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM $MEM --left$LIST $fwd --right$LIST $rev --output $dir --CPU $PROC --max_cov $KMER --pairs_together $@
else
  singularity exec $singularity \
	/usr/local/bin/util/insilico_read_normalization.pl --seqType fq --JM $MEM --single $fwd --output $dir --CPU $PROC --max_cov $KMER $@
fi
