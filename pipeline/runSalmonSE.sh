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
GC="--gcBias"
SEQ="--seqBias"
VAL="--validateMappings"
LTYPE="A"
BIND=/mnt:/mnt 
IMG=/mnt/picea/projects/singularity/salmon.simg

# USAGE
## a usage function
USAGETXT=\
"
    Usage: $0 [options] <index dir> <fq file> <out dir> <fragment length mean> <fragment length sd>
    
    Options:
              b: bind directory path
              i: path to the singularity image
              g: turn off GC bias correction (on by default)
              l: the library type, defaults to A (automatic) - see https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype 
              p: number of CPU to use (default 2)
              s: turn off sequence bias correction (on by default)
              v: turn off validateMappings (on by default), faster, less accurate
" 

# Check the tool is avail
# we use singularity
# isExec salmon

## get the options
while getopts b:i:gl:p:sv option
do
    case "$option" in
    b) BIND="$OPTARG";;
    i) IMG="$OPTARG";;
    g) GC="";;
    l) LTYPE="";;
	  p) CPU="$OPTARG";;
	  s) SEQ="";;
	  v) VAL="";;
	  \?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

## we get 2 files, 1 dir and 2 numbers as input 
if [ $# != 5 ]; then
   abort "This function has 5 arguments"
fi

if [ ! -d $1 ]; then
   abort "The first argument needs to be an existing index directory"
fi

if [ ! -f $2 ]; then
   abort "The second argument needs to be an existing fq file"
fi

if [ ! -d $3 ]; then
   abort "The third argument needs to be an existing directory"
fi

if [ $4 -lt 1 ]; then
   abort "The forth argument needs to be a positive number"
fi

if [ $5 -lt 1 ]; then
   abort "The fifth argument needs to be a positive number"
fi


bindDir=$(echo $BIND | cut -d: -f1)
if [ ! -d $bindDir ]; then
   abort "The bind dir needs to be an existing directory"
fi

if [ ! -f $IMG ]; then
  abort "The singularity image needs to be an existing file"
fi

# create an output dir based on the sample name
outdir=$3/$(basename ${2/.f*q.gz/})

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

# set the options
OPTIONS="$GC $SEQ $VAL --fldMean $4 --fldSD $5"

# run
 singularity exec --bind $BIND $IMG \
 salmon quant -p $CPU -i $1 -l $LTYPE -r $2 -o $outdir $OPTIONS 
