#! /bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 0-00:10:00
#SBATCH --mail-type=ALL

##
set -e

## 
echo Loading
export LD_LIBRARY_PATH=~delhomme/lib

##
echo Checking

## we get two dir as input
if [ $# != 4 ]; then
    echo "This function takes two directories and two files as arguments."
    echo "Usage: sbatch runCuffmerge.sh <in dir> <out dir> <gene gff3> <genome softmasked fasta>"
    exit 1
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory"
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing file"
fi

if [ ! -f $4 ]; then
    echo "The forth argument needs to be an existing file"
fi

##
echo Starting

cuffcompare -r $3 -R -C -V -s $4 -o $2 $1/*/*_transcripts.gtf > $2/cuffcompare.txt 2> $2/cuffcompare.err

##
echo Done

