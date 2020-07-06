#!/bin/bash
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH --mail-type ALL

# fail on ERROR
set -eux

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

CPU=8
OPTIONS=
IMG=/mnt/picea/projects/singularity/salmon.simg
DECOY=

# usage
USAGETXT=\
"
  $0 [options] <transcript file> <output file>

  Options:
  -d the genome file to use to build the decoy sequences
  -i the salmon image to use, defaults to salmon.simg
  -t number of threads (default 8)
  -p triggers --perfectHash
"

# process the arguments
## get the options
while getopts d:i:t:p option
do
  case "$option" in
      d) DECOY=$OPTARG
        if [ ! -f $DECOY ]; then
          abort "The genome file to extract the decoys from does not exist"
        fi
        ;;
      i) IMG=$OPTARG;;
	    t) CPU=$OPTARG;;
	    p) OPTIONS="--perfectHash $OPTIONS";;
		  \?) ## unknown flag
		  abort "Unknow option";;
  esac
done
shift `expr $OPTIND - 1`

# test
if [ "$#" -ne "2" ]; then
  abort "This script expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument should be a fasta file"
fi
tx=$1
shift

if [ ! -d $1 ]; then
  abort "The output directory does not exist, create it first"
fi
out=$1
shift

if [ $DECOY != "" ]; then
  grep "^>" $DECOY | cut -d " " -f 1 | tr -d '>' > $out/decoys.txt
  cat $tx $DECOY | gzip -c > $out/gentrome.fa.gz
  OPTIONS="-d $out/decoys.txt $OPTIONS"
  tx=$out/gentrome.fa.gz
fi

# exec
singularity exec --bind /mnt:/mnt $IMG salmon index -t $tx -i $out -p $CPU $OPTIONS
