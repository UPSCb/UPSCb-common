#!/bin/bash
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH --mail-type ALL

# TODO implement that if -d
#grep "^>" ../../../fasta/Potra01-genome.fa | cut -d " " -f 1 | tr -d '>' > decoys.txt
#cat ../../../fasta/Potra01-mRNA.fa ../../../fasta/Potra01-genome.fa | gzip -c > gentrome.fa.gz

# fail on ERROR
set -eux

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

CPU=8
OPTIONS=
IMG=/mnt/picea/projects/singularity/salmon.simg

# usage
USAGETXT=\
"
  $0 [options] <transcript file> <output file>

  Options:
  -d a text file containing the IDs of decoy sequences present in the transcript fasta file
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
          abort "The decoy file does not exist"
        fi
        OPTIONS="-d $DECOY $OPTIONS";;
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

if [ ! -d `dirname $2` ]; then
  abort "The output directory does not exist, create it first"
fi

# exec
singularity exec --bind /mnt:/mnt $IMG salmon index -t $1 -i $2 -p $CPU $OPTIONS
