#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
## SBATCH --mem=24G

## stop on error and be verbose
set -ex

# ## load the modules
module load bioinfo-tools gmap-gsnap

## usage
usage(){
echo >&2 \
"
	Usage: runGMAPIndex.sh <index dir> <index name> <fasta file>
"
	exit 1
}

## we get one dir, one token and one file as input
if [ $# != 3 ]; then
    echo "This function takes one directory, one token and one file as arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the GMAP index directory"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be a fasta file"
    usage
fi

## run GMAP
echo Indexing

## run
gmap_build -D $1 -d $2 $3

## fix permission
chmod -R g+w $1/$2

##
echo Done


