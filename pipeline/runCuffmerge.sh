#! /bin/bash -l
#SBATCH -p node -N 1
#SBATCH -t 0-02:00:00
#SBATCH --mail-type=ALL

##
set -e

## 
echo Loading
export LD_LIBRARY_PATH=~delhomme/lib

##
echo Checking

## we get two dir as input
if [ $# -lt 3 ]; then
    echo "This function takes two directories and two files as arguments, the second file being facultative"
    echo "in which case cuffmerge is run without prior annotation knowledge."
    echo "Usage: sbatch runCuffmerge.sh <in dir> <out dir> <genome fasta> [gene gff3]"
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

ext=
if [ ! -z $4 ]; then
    if [ ! -f $4 ]; then
	echo "The forth argument needs to be an existing file"
    fi
    ext="--ref-gtf $4"    
fi

##
echo Setting up

find $1 -name "transcripts.gtf" > $2/cuffmergeMANIFEST

##
echo Starting

echo cuffmerge -p 8 $ext -s $3 -o $2 $2/cuffmergeMANIFEST

##
echo Done

