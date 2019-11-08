#!/bin/bash -l
#SBATCH -p core
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00

# stop on error, be verbose and expand the commands
set -e -x

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# options
PREP=0
LEFT=
RIGHT=

## usage
USAGETXT=\
"
	Usage: runTrinityEstimeExpression.sh [options] <Trinity.fasta> <out dir>
	
	Options: -p prepare the reference (salmon index). This is mandatory before quantification
	         -f The forward file; ignored if -p is set
	         -r The reverse file; ignored if -p is set
"

## get the options
while getopts f:pr: option
do
  case "$option" in
      f) LEFT=$OPTARG;;
	    p) PREP=1;;
	    r) RIGHT=$OPTARG;;
		\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

# Check
if [ $# -ne 2 ]; then
    echo "This function needs 2 arguments"
    usage
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the trinity fasta filepath"
fi

if [ ! -d $2 ]; then
    abort "The second argument (output dir) needs to be an existing directory"
fi

# go in the out directory
cd $2

# run
if [ $PREP -ne 0 ]; then
  singularity exec --bind /mnt:/mnt  /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
 /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl \
 --transcripts $1 --est_method salmon --trinity_mode --prep_reference
else
  if [ ! -f $LEFT ]; then
    abort "For quantification, you need to provide an existing forward file; set the -f option to the correct file path"
  fi

  if [ ! -f $RIGHT ]; then
    abort "For quantification, you need to provide an existing reverse file; set the -r option to the correct file path"
  fi

  singularity exec --bind /mnt:/mnt  /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
 /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl \
 --transcripts $1 --est_method salmon --trinity_mode --seqType fq \
 --left $LEFT --right $RIGHT --output_dir $2
fi

