#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools
module load samtools/0.1.19

# usage 
usage(){
echo >&2 \
"
	Usage: $0 <array.tsv> <in.bam> <out dir>
        
        Purpose: This function reheader the sequence name within bam files.

        Note: The array file should contain 2 columns, one with the old ID and the second one with the new ID separated by space or tabs
              The output name remains unchanged in the output directory
"
	exit 1
}

## we get two files and one dir as input
if [ $# != 3 ]; then
    echo "This function takes two files and one dir as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing file mapping old to new sequence IDs"
    usage
fi

if [ ! -f $2 ]; then
    echo "The first argument needs to be an existing BAM file"
    usage
fi

if [ ! -d $3 ]; then
    echo "The third argument needs to be an existing directory"
    usage
fi

if [ $3 == `dirname $2` ]; then
    echo "We do not inplace edit, indir must be different from outdir"
    exit 1
fi

ARRAY=$1
## check if file is zipped
fdat=$(file $ARRAY)
if [[ $fdat =~ gzip ]]; then
    ## if yes, FIFO it with a named pipe
    CLEAN_FIFO=true
    tmpf=/tmp/tempfifo${RANDOM}.conversion
    mkfifo $tmpf
    zcat $ARRAY > $tmpf &
    ARRAY=$tmpf
fi
## reheader 
samtools view -H $2 | awk 'BEGIN {OFS="\t"} FNR==NR{a[$1]=$2;next}{if($1=="@SQ"){sub(/SN:/,"",$2);print $1, "SN:"a[$2],$3}else{print $0}}' $ARRAY - | samtools reheader - $2 > $3/`basename ${2}`
samtools index $3/`basename ${2}`
if [ ! -z CLEAN_FIFO ]; then
    rm $tmpf
fi
