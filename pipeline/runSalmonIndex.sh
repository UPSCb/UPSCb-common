#!/bin/bash
#SBATCH -p node -n 8
#SBATCH -t 12:00:00
#SBATCH --mail-type ALL

# fail on ERROR
set -eux

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

CPU=8
OPTIONS=
DECOY=

# usage
USAGETXT=\
"
  $0 [options] <singularity image> <transcript file> <output dir>

  Options:
  -d the genome file to use to build the decoy sequences
  -k the kmer size (default to salmon's default: 31)
  -t number of threads (default 8)
  -p triggers --perfectHash
"

# process the arguments
## get the options
while getopts d:k:t:p option
do
  case "$option" in
      d) DECOY=$OPTARG
        if [ ! -f $DECOY ]; then
          abort "The genome file to extract the decoys from does not exist"
        fi
        ;;
      k) OPTIONS="-k $OPTARG $OPTIONS";;
	    t) CPU=$OPTARG;;
	    p) OPTIONS="--perfectHash $OPTIONS";;
		  \?) ## unknown flag
		  abort "Unknow option";;
  esac
done
shift `expr $OPTIND - 1`

# test
if [ "$#" -ne "3" ]; then
  abort "This script expects 3 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument should be a singularity file"
fi
IMG=$1
shift

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

if [ ! -f $1 ]; then
  abort "The second argument should be a fasta file"
fi
tx=$1
shift

if [ ! -d $1 ]; then
  abort "The output directory does not exist, create it first"
fi
out=$1
shift

cd $out
if [ ! -z $DECOY ]; then
  gunzip -c $DECOY | grep "^>" | cut -d " " -f 1 | tr -d '>' > decoys.txt
  cat $tx $DECOY > gentrome.fa.gz
  OPTIONS="-d decoys.txt $OPTIONS"
  tx=gentrome.fa.gz
fi

# exec
singularity exec $IMG salmon index -t $tx -i . -p $CPU $OPTIONS
