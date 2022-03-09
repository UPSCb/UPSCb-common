#!/bin/bash -l
#SBATCH -p core -n 2
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=12GB

# be verbose and fail on errr
set -ex

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# defaults
CPU=2
ECLASS="--dumpFeatures"
STYPE="--chromiumV3"
LTYPE="ISR"
BIND=/mnt:/mnt 
IMG=/mnt/picea/projects/singularity/salmon.simg

# USAGE
## a usage function
USAGETXT=\
"
    Usage: $0 [options] <index dir> <fwd fq file> <rev fq file> <out dir> <transcript_map>
    
    Options:
              b: bind directory path
              e: dump equivalence class - turn on --dumpFeatures
              i: path to the singularity image
              l: the library type, defaults to ISR- see https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype 
              p: number of CPU to use (default 2)
" 

# Check the tool is avail
# we use singularity
# isExec salmon

## get the options
while getopts b:e:i:l:p: option
do
    case "$option" in
    b) BIND="$OPTARG";;
    e) ECLASS="--dumpFeatures";;
    i) IMG="$OPTARG";;
    l) LTYPE="$OPTARG";;
	  p) CPU="$OPTARG";;
	  \?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

# set the options
OPTIONS="$STYPE $ECLASS"

## we get 4 files and 1 dir as input 
if [ $# != 5 ]; then
   abort "This function has 5 arguments"
fi

if [ ! -d $1 ]; then
   abort "The first argument needs to be an existing index directory"
fi

if [ ! -f $2 ]; then
   abort "The second argument needs to be an existing fq file"
fi

if [ ! -f $3 ]; then
   abort "The third argument needs to be an existing fq file"
fi

if [ ! -d $4 ]; then
   abort "The fourth argument needs to be an existing directory"
fi

if [ ! -f $5 ]; then
   abort "The fifth argument needs to be an existing transcript map file"
fi

bindDir=$(echo $BIND | cut -d: -f1)
if [ ! -d $bindDir ]; then
   abort "The bind dir needs to be an existing directory"
fi

if [ ! -f $IMG ]; then
  abort "The singularity image needs to be an existing file"
fi

# create an output dir based on the sample name
outdir=$4/$(basename ${2/_1.f*q.gz/})

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

# run
 singularity exec --bind $BIND $IMG \
 salmon alevin -p $CPU -i $1 -l $LTYPE -1 $2 -2 $3 -o $outdir $OPTIONS --tgMap $5