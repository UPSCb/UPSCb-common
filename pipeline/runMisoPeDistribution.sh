#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 0-03:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

## samtools
module load bioinfo-tools
module load samtools/0.1.19

## check for the utility
MISO=`which pe_utils.py`
if [ $? != 0 ]; then
    echo "please install MISO before running this script or add it to your PATH"
    exit 1
fi

if [ ! -f $MISO -a ! -x $MISO ]; then
    echo "your MISO pe_utils.py does not appear to be an executable file"
    exit 1
fi

MISOPLOT=`which plot.py`
if [ $? != 0 ]; then
    echo "please install MISO before running this script or add it to your PATH"
    exit 1
fi

if [ ! -f $MISOPLOT -a ! -x $MISOPLOT ]; then
    echo "your MISO plot.py does not appear to be an executable file"
    exit 1
fi

## usage
usage(){
echo >&2 \
"
	Usage: runMisoPeDistribution.sh <out dir> <bam file> <gff file>
"
exit 1
}

## check the args
if [ ! -d $1 ]; then
    echo "The output directory: $1  does not exist"
    usage
fi

if [ ! -f $2 ]; then
    echo "The input bam file: $2  does not exist"
    usage
fi

if [ ! -f $3 ]; then
    echo "The gff file: $3  does not exist"
    usage
fi

## run
$MISO --compute-insert-len $2 $3 --output-dir $1

## plot
echo "[data]
bam_prefix = `dirname $2` 
miso_prefix = $1
bam_files = [ '`basename $2`']
miso_files = [ '`basename $2`.insert_len']" > $2-sashimi.settings

$MISOPLOT --plot-insert-len $1/`basename $2`.insert_len $2-sashimi.settings --output-dir $1

