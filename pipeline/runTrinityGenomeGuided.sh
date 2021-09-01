#!/bin/bash -l
#SBATCH -p core
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 16
#SBATCH --mem=128G
## time too for large files
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

## stop on error, be verbose and expand the commands
set -eux

## load necessary modules
#module load bioinfo-tools
#module load bowtie
#module load samtools
#module load trinity

## check the options if any
PROC=20
INTRON=11000
MEM="128G"
STSPEC=
LONGREADS=
NORM=

# source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
USAGETXT=\
"
	Usage: runTrinityGenomeGuided.sh [options] <out dir> <BAM file>
	
	Options:
	              -i max intron size
                -l long reads (pac bio or other) as a fasta file
                -m mem requirement (in GB, default 128)
                -n no digital read normalisation
		            -p number of threads to use (default 16)
                -s strand specific (Illumina dUTP protocol)
"

## get the options
while getopts l:m:np:si: option
do
    case "$option" in
      i) INTRON=$OPTARG;; 
      l) LONGREADS="$OPTARG";;
      n) NORM="--no_normalize_reads";;
      m) MEM=$OPTARG;;
      p) PROC=$OPTARG;;
      s) STSPEC="--SS_lib_type RF";;
	 	\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

# sanity
## we get one dir and one file
if [ $# != 2 ]; then
    abort "This function needs 2 arguments"
fi

if [ ! -d $1 ]; then
    abort "The first argument (output dir) needs to be an existing directory"
fi

if [ ! -f $2 ]; then
    echo "The second argument (alignment BAM) needs to be an existing file"
    usage
fi

## check args
if [ ! -z $LONGREADS ]; then
	if [ ! -f $LONGREADS ]; then
		abort "-l should point to an existing fasta file"
	else
		LONGREADS="--long_reads $LONGREADS"
	fi	
fi


## run trinity
echo Assembling
#$TRINITY_HOME/Trinity --genome_guided_bam $5 --genome_guided_max_intron $INTRON --CPU $GGCPU --seqType fq --max_memory ${JMEM}G  --left $2 --right $3 --CPU $PROC --output $1 --min_kmer_cov $KMER $STSPEC
singularity exec --bind /mnt:/mnt -e /mnt/picea/projects/singularity/trinityrnaseq.v2.13.1.simg \
Trinity --genome_guided_bam $2 --genome_guided_max_intron $INTRON \
--max_memory $MEM --CPU $PROC --output $1 $LONGREADS $NORM $STSPEC

##
echo Done

