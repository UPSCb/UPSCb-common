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
BC=
OUTCL=4

## usage
USAGETXT=\
"
	Usage: runJellyfishBc.sh [options] <in.bam> <out.dir>

	Options:
	  -b the path to a bc file
	  -C turn off canonical mode (on by default)
	  -l set --out-counter-len (default 4, can be set to 2 or 1 to reduce mem usage)
	  -m kmer size (default $KMERSIZE)
	  -s hash size (defaut $HASHSIZE)
	  -t threads (default $CPUS)
"

while getopts b:Cl:m:s:t: option
do
        case "$option" in
        b) BC=$OPTARG;;
        C) CANON="";;
        l) OUTCL=$OPTARG;;
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

[[ ! -z $BC ]] && [[ ! -f $BC ]] && abort "The -b argument should be an existing file"

[[ ! -z $BC ]] && BC="--bc $BC"

[[ $OUTCL -ne 1 ]] || [[ $OUTCL -ne 2 ]] || [[ $OUTCL -ne 4 ]] && abort "-l must be a value among 1,2 or 4"

isExec jellyfish
isExec samtools

jellyfish count -m $KMERSIZE -s $HASHSIZE $CANON \
-o $2/$(basename ${1/.bam/.jf}) -t $CPUS $BC --out-counter-len $OUTCL <(samtools view $1 | awk '{print ">"$1"\n"$10}')

