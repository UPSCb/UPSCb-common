#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 16
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
	Usage: $0 [option] <in.bam>
	Options:
                -p define the number of threads to use (16)
                -n sort by name instead of coordinates
                -i inplace sort - i.e. keep the file name unchanged
        Note: If sorting by coordinates, the extension is _byCoord.bam
              and for sorting by ID, the extenstion is _byReadID.bam
"
	exit 1
}

# define options
SORT=
EXT="_byCoord"
INPLACE=0
CPU=16
## get the options
while getopts p:ni option
do
        case "$option" in
	    i) INPLACE=1;;
	    n) SORT=" -n "
		EXT="_byReadID"
		;;
	    p) CPU=$OPTARG;;
	    \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get one file as input
if [ $# != 1 ]; then
    echo "This function takes one file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing bam file"
fi

## define the output file
new=`dirname $1`/`basename ${1//.bam/}`${EXT}

## get the coverage table
samtools sort -@ $CPU $SORT $1 $new

## if inplace
if [ $INPLACE == 1 ]; then
    mv $new.bam $1
fi

