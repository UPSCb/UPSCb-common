#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

## module load
module load bioinfo-tools
module load samtools/0.1.19
module load python/2.7.1

## usage
usage(){
echo >&2 \
"
	Usage: runMiso.sh <output dir> <index dir> <bam file> 

	Options:
                -l read length (default to 101)
                -m mean fragment size (default to 250)
                -s fragment size sd (default to 25)
"
exit 1
}

## get the options
LEN=101
MEAN=250
SD=25
while getopts l:m:s: option
do
        case "$option" in
	    l) LEN=$OPTARG;;
	    m) MEAN=$OPTARG;;
	    s) SD=$OPTARG;;
	    \?) ## unknown flag
	    usage;;
        esac
done
shift `expr $OPTIND - 1`

## check for the utilities
MISO=`which run_events_analysis.py`
if [ $? != 0 ]; then
    echo "please install MISO before running this script or add it to your PATH"
    usage
fi

if [ ! -f $MISO -a ! -x $MISO ]; then
    echo "your MISO run_events_analysis.py does not appear to be an executable file"
    usage
fi

MISOPY=`which run_miso.py`
if [ $? != 0 ]; then
    echo "please install MISO before running this script or add it to your PATH"
    usage
fi

if [ ! -f $MISOPY -a ! -x $MISOPY ]; then
    echo "your MISO run_miso.py does not appear to be an executable file"
    usage
fi

MISOZIP=`which miso_zip.py`
if [ $? != 0 ]; then
    echo "please install MISO before running this script or add it to your PATH"
    usage
fi

if [ ! -f $MISOZIP -a ! -x $MISOZIP ]; then
    echo "your MISO miso_zip.py does not appear to be an executable file"
    usage
fi

## check the args
if [ ! -d $1 ]; then
    echo "The output directory: $1  does not exist"
    usage
fi

if [ ! -d $2 ]; then
    echo "The index directory: $2  does not exist"
    usage
fi

if [ ! -f $3 ]; then
    echo "The input bam file: $3  does not exist"
    usage
fi

## run MISO
$MISO --compute-genes-psi $2 $3 --output-dir $1 --prefilter --read-len $LEN --paired-end $MEAN $SD -p 8

## summarize
$MISOPY --summarize-samples $1 $1/summary

## compress
## $MISOZIP --compress $1.misozip $1

## remove some data
## rm $1/*.miso

