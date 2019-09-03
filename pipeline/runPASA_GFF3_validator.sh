#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL

## stop on error, be verbose
set -ex

## load modules
module load bioinfo-tools pasa

## a usage function
usage(){
    echo >&2 "Usage: $0 <gff file>" 
    exit 1
}

## check the number of arguments
if [ $# != 1 ]; then
    echo "This function take one argument, the gff file path"
    usage
fi

## check that the file exists
if [ ! -f $1 ]; then
    echo "The argument shoud be a valid file path to a gff file"
    usage
fi

## run it
perl $PASAHOME/misc_utilities/pasa_gff3_validator.pl $1
