#!/bin/bash -l
#SBATCH -p core -n 12
#SBACTH --mail-type=FAIL
#SBATCH -t 2-00:00:00

set -eux

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# vars
OPTARG=
OPTIND=
CANON="-C"
OPTIONS=
KMERSIZE=25
HASHSIZE=100G
CPUS=12

## usage
USAGETXT=\
"
	Usage: runJellyfishBc.sh <in.bam> <out.dir>
	
	Options:
	  -C turn off canonical mode (on by default)
	  -m kmer size (default $KMERSIZE)
	  -s hash size (defaut $HASHSIZE)
	  -t threads (default $CPUS)
"

while getopts Cm:s:t: option
do
        case "$option" in
        C) CANON="";;
	      m) KMERSIZE=$OPTARG;;
        s) HASHSIZE=$OPTARG;;
        t) CPUS=$OPTARG;;
        \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

# sanity
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file"

[[ ! -d $2 ]] && abort "The second argument needs to be an existing directory"

isExec jellyfish
isExec samtools

cd $2

jellyfish bc -m $KMERSIZE -s $HASHSIZE $CANON -o $(basename ${1/.bam/.bc}) -t $CPUS <(samtools view $1 | awk '{print ">"$1"\n"$10}')

