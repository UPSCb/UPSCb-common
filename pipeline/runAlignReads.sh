#! /bin/bash -l
#SBATCH -p node -n 16
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## load modules
module load bioinfo-tools
module load bowtie/0.12.9
module load samtools

## vars
PROC=16

## usage
usage(){
echo >&2 \
"
	Usage: runAlignReads.sh <out dir> <tmp dir> <trinity fa> <left fq.gz> <right fq.gz>

        Note:
             You need to set the TRINITY_RNASEQ_ROOT env. variable to your trinity directory.
"
	exit 1
}

## we get two dir and three files as input
if [ $# != 5 ]; then
    echo "This function requires five  arguments."
    usage
fi

## check the out dir
if [ ! -d $1 ]; then
    echo "The first argument should be the (existing) output directory"
    usage
fi

## check the tmp dir
if [ ! -d $2 ]; then
    echo "The second argument should be the (existing) tmp directory"
    usage
fi

## check the files
if [ ! -f $3 ]; then
    echo "The third argument should be the trinity fasta file"
    usage
fi

if [ ! -f $4 ]; then
    echo "The forth argument should be the forward fq.gz file"
    usage
fi


if [ ! -f $5 ]; then
    echo "The fifth argument should be the reverse fq.gz file"
    usage
fi

## check for TRINITY_RNASEQ_ROOT
if [ -z $TRINITY_RNASEQ_ROOT ]; then
    echo "The TRINITY_RNASEQ_ROOT env. var. should be set to your Trinity install dir."
    usage
fi

## decompress
printf "$4\0$5" | xargs -0 -P 2 -I {} bash -c 'gunzip -c $0 > $1/`basename ${0//.gz/}`' {} $2

## create the out dir and cd there. The script does print files in the current dir...
out=$1/`basename ${4//_1.f*q.gz/}`

if [ ! -d $out ]; then
    mkdir $out
fi

cd $out

## realign
$TRINITY_RNASEQ_ROOT/util/align_and_estimate_abundance.pl --transcripts $3 --left $2/`basename ${4//.gz/}` --right $2/`basename ${5//.gz/}` --seqType fq --aln_method bowtie --est_method RSEM --output_dir $out --thread_count $PROC

