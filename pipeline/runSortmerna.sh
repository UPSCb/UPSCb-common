#!/bin/bash
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
## time too for large files
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

# sanity check
set -eux

# check the options if any
#KEEP=1
UNPAIRED=0
PROC=20
INX=
SEEDS=2
PASSES=
SAM=
DB=/mnt/picea/storage/reference/rRNA/sortmerna/v4.2

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# alias the exec kogia container
sortmerna=${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/kogia/scripts/sortmerna

# usage
export USAGETXT="
	Usage: runSortmerna.sh [option] <out dir> <forward fastq.gz> <reverse fastq.gz>

	Options:
                -a report alignments (default is not to report)
                -d the database directory
                -i index dir
                -p number of threads to be used (default $PROC)
                -P number of passes (default to sortmerna defaults)
		            -s number of seeds (default 2)
                -u single end data (in that case only the forward fastq is needed)

  Note:
    All the fasta files found in the <db dir> will be used. If -i is provided this will be 
    passed as an argument to SortMeRNA in addition.
"
# get the options
while getopts ad:i:p:P:s:u option
do
  case "$option" in
        a) SAM="--sam --num_alignments 1";;
        d) DB=$OPTARG;;
        i) INX="--idx $OPTARG";;
  	    p) PROC=$OPTARG;;
	      P) PASSES="--passes $OPTARG";;
	      s) SEEDS=$OPTARG;;
	      u) UNPAIRED=1;;
		    \?) ## unknown flag
		    abort;;
  esac
done
shift `expr $OPTIND - 1`

# we get two dir and two files as input
if [ $UNPAIRED == 0 ]; then
    if [ $# != 3 ]; then
	    abort "This function takes one directory and two files as arguments"
    fi
    if [ ! -f $3 ]; then
      abort "The fourth argument needs to be an existing file"
    fi
else
    if [ $# != 2 ]; then
	    abort "This function takes one directory and one file as argument"
    fi
fi

if [ ! -d $1 ]; then
    abort "The first argument needs to be an existing directory"
fi

if [ ! -d $2 ]; then
    abort "The second argument needs to be an existing directory"
fi

if [ ! -f $3 ]; then
    abort "The third argument needs to be an existing file"
fi

# PE
if [ $UNPAIRED == 0 ]; then
    fo=$(basename ${3//_[1,2].f*q.gz/_sortmerna})
else
    fo=$(basename ${3//.f*q.gz/_sortmerna})
fi

# DBs
dbs=$(find $DB -maxdepth 1 -mindepth 1 -name "*.fasta" -type f -printf " --ref %p")

# run
if [ $UNPAIRED == 0 ]; then
    sortmerna $dbs --workdir $1/$fo \
    --reads $3 --reads $4 $IDX --threads $PROC \
    --fastx --paired_in --out2 --other $1/$fo $INX \
    --aligned $1/${fo}_rRNA $SAM $PASSES --num_seeds $SEEDS 
else
    sortmerna $dbs --reads $3 $IDX --threads $PROC \
    --fastx --other $1/$fo --workdir $1 \
    --aligned $1/${fo}_rRNA $SAM $PASSES --num_seeds $SEEDS 
fi

## rm the out
rm -rf $1/$fo

## compress the output files
find $1 -name "${fo}*.fastq" -print0 | xargs -0 -I {} -P $PROC gzip -f {}
