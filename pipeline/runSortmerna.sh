#!/bin/bash
#SBATCH -p main
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=6GB

# stop on error
set -eu

# be verbose and extend the commands
#set -x

# check the options if any
UNPAIRED=0
PROC=20
DBS=
SEEDS=
PASSES=
SAM=
ILV=0

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
export USAGETXT="
	Usage: runSortmerna.sh [option] <singularity> <out dir> <ref fasta> <inx dir> <forward fastq.gz> <reverse fastq.gz>

	Options:
                -a report alignments (default is off)
                -p number of threads to be used (default $PROC)
                -P number of passes (default to sortmerna defaults)
		            -s number of seeds (default 2)
                -u single end data (in that case only the forward fastq is needed)
            		-i interleaved data (in that case only one fastq is needed)
"
# get the options
while getopts aip:P:s:u option
do
        case "$option" in
        a) SAM="--sam --num_alignments 1";;
	      i) ILV=1
	        UNPAIRED=1;;
  	    p) PROC=$OPTARG;;
	      P) PASSES="--passes $OPTARG";;
	      s) SEEDS="--num_seeds $OPTARG";;
	      u) UNPAIRED=1;;
		    \?) ## unknown flag
		    abort;;
        esac
done
shift `expr $OPTIND - 1`

# we get two dir and two files as input
if [ $UNPAIRED == 0 ]; then
    [[ $# != 6 ]] && abort "This function takes two directories and four files as arguments"
else
    [[ $# != 5 ]] &&abort "This function takes two directories and three files as argument"
fi

[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity sortmerna container file"

[[ ! -d $2 ]] && abort "The second argument needs to be an existing output directory"

[[ ! -f $3 ]] && abort "The third argument needs to be an existing rRNA fasta file"

[[ ! -d $4 ]] && abort "The fourth argument needs to be an existing rRNA index directory"

[[ ! -f $5 ]] && abort "The fifth argument needs to be an existing fastq.gz file"

if [ $UNPAIRED == 0 ]; then
    [[ ! -f $6 ]] && abort "The sixth argument needs to be an existing fastq.gz file"
fi

# sample prefix
fnam=$(basename ${5/.f*q.gz/})
out=$2/$fnam
[[ ! -d $out ]] && mkdir -p $out

# options
OPTIONS="$PASSES $SEEDS --threads $PROC"
if [ ! -z $SAM ]; then
  OPTIONS="$OPTIONS $SAM --aligned $out/${fnam}_rRNA"
fi

## rm the tmp
if [ $UNPAIRED == 0 ]; then
    singularity exec $1 sortmerna --ref $3 --workdir $out --idx-dir $4 --reads $5 --reads $6 --other $out/${fnam}_sortmerna --fastx --paired_in --out2 $OPTIONS
else
  if [ $ILV -eq 1 ]; then
    singularity exec $1 sortmerna --ref $3 --workdir $out --idx-dir $4 --reads $5 --other $out/${fnam}_sortmerna --paired_in --fastx --out2 $OPTIONS
  else
    singularity exec $1  sortmerna --ref $3 --workdir $out --idx-dir $4 --reads $5 --other $out/${fnam}_sortmerna --fastx $OPTIONS
  fi
fi
