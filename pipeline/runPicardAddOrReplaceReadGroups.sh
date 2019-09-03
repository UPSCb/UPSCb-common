#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mem=32GB

## stop on error
set -e

## modules
module load bioinfo-tools Picard-tools samtools

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 <bam file> <out dir> <read groups>

    Note:
         Read Groups should be provided as ID=value pairs space delimited
" 
    exit 1
}

## we get one file as input
if [ $# -le 3 ]; then
    echo "This function takes at least a bam file, an out dir and one read group specification as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing bam file"
    usage
fi
in=$1
shift

if [ ! -d $1 ]; then
    echo "The second argument needs to be an existing directory"
    usage
fi
out=$1
shift

# create the outfile
out=$out/`basename $in`

## clean sam
java -jar $PICARD_TOOLS_DIR/picard.jar CleanSam I=$in O=$out.clean

## add the RG
java -jar $PICARD_TOOLS_DIR/picard.jar AddOrReplaceReadGroups I=$out.clean O=$out $@

# index the bam file
samtools index $out

## clean up
rm $out.clean
