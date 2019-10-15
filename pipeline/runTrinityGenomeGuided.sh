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
set -e -x

## load necessary modules
module load bioinfo-tools
#module load bowtie
#module load samtools
module load trinity

## check the options if any
PROC=16
INTRON=11000
KMER=1
MEM="128G"
STSPEC=

## usage
usage(){
echo >&2 \
"
	Usage: runTrinityGenomeGuided.sh [options] <out dir> <left fq(s)> <right fq(s)> <BAM file>
	
	Options:
                -k min kmer cov (default 1. Change to 2 w/o digital normalization) 
                -m mem requirement (in GB, default 128)
		            -p number of threads to use (default 16)
                -s strand specific (Illumina dUTP protocol)
		-i max intron size

        Note:
             You need to set the TRINITY_HOME env. variable to your trinity directory.
"
	exit 1
}

## get the options
while getopts k:m:p:si: option
do
    case "$option" in
      k) KMER=$OPTARG;;
      m) MEM=$OPTARG;;
      p) PROC=$OPTARG;;
      s) STSPEC="--SS_lib_type RF";;
      i) INTRON=$OPTARG;; 
	 	\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function needs 4 arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument (output dir) needs to be an existing directory"
    usage
fi

## proceed the other args
if [ ! -f $2 ]; then
    echo "The second argument (forward reads) needs to be an existing file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument (reverse reads) needs to be an existing file"
    usage
fi

if [ ! -f $4 ]; then
    echo "The forth argument (alignment BAM) needs to be an existing file"
    usage
fi

## run trinity
echo Assembling
#$TRINITY_HOME/Trinity --genome_guided_bam $5 --genome_guided_max_intron $INTRON --CPU $GGCPU --seqType fq --max_memory ${JMEM}G  --left $2 --right $3 --CPU $PROC --output $1 --min_kmer_cov $KMER $STSPEC
$TRINITY_HOME/Trinity --genome_guided_bam $4 --genome_guided_max_intron $INTRON \
--seqType fq --max_memory $MEM --left $2 --right $3 --CPU $PROC --output $1 \
--min_kmer_cov $KMER $STSPEC

##
echo Done

