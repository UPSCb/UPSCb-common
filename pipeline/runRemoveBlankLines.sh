#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## usage
usage(){
echo >&2 \
"
	Usage: $0 <gz file>
"
	exit 1
}

## we get file as input
if [ $# != 1 ]; then
    echo "This function takes one file as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing gzipped file"
    usage
fi

# procced
zcat $1 | sed '/^$/d' | gzip > $1.tmp
mv $1.tmp $1
