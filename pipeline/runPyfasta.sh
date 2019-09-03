#! /bin/bash -l
#SBATCH -p node -n 2
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

## we need slightly more than 3GB of RAM :-\

## stop on error
set -e

if [ $# != 3 ]; then
    echo "The argument should be the fasta full filename, the output directory and the number of chunks"
    exit 1
fi

if [ ! -f $1 ]; then
    echo "The  fasta filename you provided does not exist"
    exit 1
fi

if [ ! -d $2 ]; then
    echo "The  directory name you provided does not exist"
    exit 1
fi

if [ -z $3 ]; then
    echo "The second argument should be an integer value describing the final number of chunks expected"
    exit 1
fi

if [ "$3" -lt "1" ]; then
    echo "The second argument should be an integer value larger than 1 describing the final number of chunks expected"
    exit 1
fi

## cd in the out dir
cd $2

## copy if it does not exist
if [ ! -f  $2/`basename $1` ]; then
cp $1 $2
fi

## create the chunks
pyfasta split -n $3 $2/`basename $1`

## rm the file
rm $2/`basename $1`


