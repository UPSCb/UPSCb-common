#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH -p main -n 1
#SBATCH -t 04:00:00
#SBATCH --mem=8G

## failsafe
set -eu

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## Vars
GENOME=12000000
NAME=
CONTROL=1
MODE="-f BAM"
CUTOFF="--cutoff-analysis"

# usage
USAGETXT=\
"
	Usage: runMacs2.sh [options] <singularity container> <treatment file> <control file> <out dir> -- [additional MACS arguments]

  Options:
    -a: cutoff analysis is on by default, set to a numeric value to deactivate and set your own MACS2 -p threshold
    -c: no control file
    -g: set the genome size; default to 12,000,000
    -n: set the output name
    -p: paired-end mode
    --: a special arguments, arguments passed afterwards are forwarded as-is to MACS2
"

# get the options
while getopts a:cg:n:p option
do
        case "$option" in
        a) CUTOFF="-p $OPTARG";;
        c) CONTROL=0;;
        g) GENOME=$OPTARG;;
        n) NAME="-n $OPTARG";;
      	p) MODE="-f BAMPE";;
		    \?) # unknown flag
		      usage;;
        esac
done
shift `expr $OPTIND - 1`


# sanity
# args
[[ $CONTROL -eq 0 ]] && [[ $@ -lt 3 ]] && abort "This function takes two files and a dir as arguments"

[[ $CONTROL -eq 1 ]] && [[ $@ -lt 4 ]] && abort "This function takes three files and a dir as arguments"

# singularity
macs2=$1
shift
[[ ! -f $macs2 ]] && abort "The first argument needs to be MACS2 singularity container"

# enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# input
treatment=$1
shift
[[ ! -f $(readlink -f $treatment) ]] && abort "The second argument needs to be the treatment BAM file"

# control
control=
if [ $CONTROL -eq 1 ]; then
  [[ ! -f $1 ]] && abort "The third argument needs to be the control bam file"
  control="-c $1"
  shift
fi

# out
out=$1
shift
[[ ! -d $out ]] && abort "The last argument needs to be the output directory"

# more args? drop --
[[ $# != 0 ]] && shift

# name
if [ -z $NAME ]; then
	if [ $CONTROL -eq 1 ]; then
	  NAME="-n $(basename ${treatment//.bam/})-$(basename ${control//.bam/})"
	else
	  NAME="-n $(basename ${treatment//.bam/})"
	fi
fi

# run
singularity exec $macs2 macs2 callpeak -t $treatment $control -g $GENOME \
--keep-dup auto $CUTOFF --outdir $out $NAME -B --SPMR $MODE  $@
