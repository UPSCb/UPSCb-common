#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e

if [ $# == 0 ]; then
    echo "This function takes one file as argument"
    exit 1
fi

if [ ! -f $1 ]; then
    echo "The provided file: $1 does not exist"
    exit 1
fi

gzip $1

