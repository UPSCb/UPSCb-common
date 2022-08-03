#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH -p rbx
#SBATCH -n 24
#SBATCH -t 12:00:00

# failsafe
set -eu

# load helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
    Usage: runBwamem.sh [options] <bwa singularity container> <samtools singularity container> <Forward read> <Reverse read> <Index Fasta> <Output dir> [--] [additional BWA arguments] 
    Options:
    -c           ----- csi index for large genomes
    -s           ----- single end alignment
    -t <INT>     ----- Number of mapping and sorting threads
    -- is a special argument that stop the command line scanning for the script  options.
 
    NOTE: You need to reserve memory for -t <INT> * 768M in SLURM
"

# vars
CSI=
SINGLE=0
THREADS=24
ARGSN=6

# options
while getopts cst: opt
do
  case "$opt" in
    c) CSI="-c";;
    s) SINGLE=1
    ARGSN=$(expr $ARGSN - $SINGLE);;
	  t) THREADS=$OPTARG;;
	  \?)# unknown flag
	        usage;;
  esac
done

shift `expr $OPTIND - 1`

# sanity
[[ $# -lt $ARGSN ]] && abort "The script expects $ARGSN arguments."

# get the containers
bwa=$1
shift
[[ ! -f $bwa ]] && abort "The first argument needs to be an existing BWA singularity container"

samtools=$1
shift
[[ ! -f $samtools ]] && abort "The first argument needs to be an existing samtools singularity container"

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# input
INPUT=$1
NAM=$(basename ${INPUT/.f*q.gz/})
shift

[[ ! -f $INPUT ]] && abort "Forward FASTQ not found"

if [ $SINGLE -eq 0 ]; then
  [[ ! -f $1 ]] && abort "Reverse FASTQ not found"
  INPUT="$INPUT $1"
  shift
fi

INDEX=$1
shift
[[ ! -f $INDEX ]] && abort "Index FASTA not found"

# output
OUT=$1
shift
[[ ! -d $OUT ]] && abort "Output directory not found"

## do we have more arguments? drop the --
[[ $# != 0 ]] && shift

# run
singularity exec $bwa bwa mem $@ -t $THREADS $INDEX $INPUT | \
singularity exec $samtools samtools view -bT $INDEX - | \
singularity exec $samtools samtools sort -@ $THREADS -o ${OUT}/${NAM}.sorted.bam

singularity exec $samtools samtools index $CSI ${OUT}/${NAM}.sorted.bam
