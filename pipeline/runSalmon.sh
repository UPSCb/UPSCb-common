#!/bin/bash -l
#SBATCH -p core -n 2
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=12GB

# be verbose and fail on errr
set -ex

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# defaults
CPU=2
ECLASS="--dumpEq --numGibbsSamples 100"
GC="--gcBias"
LTYPE="A"
POS="--posBias"
SEQ="--seqBias"
CLEAN=0

# USAGE
## a usage function
USAGETXT=\
"
    Usage: $0 [options] <singularity container> <index dir> <fwd fq file> <rev fq file> <out dir>

    Options:
              C: remove the input files upon completion of the quantification (DANGER!)
              e: turn off dump equivalence class and perform gibbs sampling (on by default)
              g: turn off GC bias correction (on by default)
              l: the library type, defaults to A (automatic) - see https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype 
              p: turn off position bias correction (on by default)
              s: turn off sequence bias correction (on by default)
              t: number of CPU to use (default 2)
"

## get the options
while getopts Cegl:pst: option
do
    case "$option" in
    C) CLEAN=1;;
    e) ECLASS="";;
    g) GC="";;
    l) LTYPE="";;
	  p) POS="";;
	  s) SEQ="";;
	  t) CPU="$OPTARG";;
	  \?) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

# set the options
OPTIONS="$GC $SEQ $POS $ECLASS"

# checks
[[ $# != 5 ]] && abort "This function expects 5 arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be the salmon singularity container file"

## enforce singularity
[[ -z $SINGULARITY_BINDPATH ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -d $2 ]] && abort "The second argument needs to be an existing index directory"

[[ ! -f $3 ]] && abort "The third argument needs to be an existing fq file"

[[ ! -f $4 ]] && abort "The fourth argument needs to be an existing fq file"

[[ ! -d $5 ]] && abort "The fifth argument needs to be an existing directory"

# create an output dir based on the sample name
outdir=$5/$(basename ${3/_1.f*q.gz/})
[[ ! -d $outdir ]] && mkdir -p $outdir

# run
singularity exec $1 salmon quant -p $CPU -i $2 -l $LTYPE -1 $3 -2 $4 -o $outdir $OPTIONS 

# clean
[[ $CLEAN -eq 1 ]] && [[ -f $outdir/quant.sf ]] && rm $3 $4
