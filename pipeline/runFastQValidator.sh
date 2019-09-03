#!/bin/bash -l

#SBATCH -p core -n 1
#SBATCH -t 0-01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16GB

## load the module if it exists
module load bioinfo-tools && module load fastQvalidator || {
  if ! hash fastQValidator 2>/dev/null; then
    echo "fastQValidator was not found in your path" 1>&2
    exit 1
  fi
}

usage() {
  echo "usage: `basename $0` <fastq>

Run fastQValidator on a FASTQ file. Prints output on stdout and
exits with a non-zero exit status if the input file does not
conform to the standard.

ARGUMENTS:
    fastq   a FASTQ file, can be gzipped

NOTES:
    fastQValidator must lie in your PATH" 1>&2
  exit 1
}

## stop on error
set -e

## check
if [ $# != 1 ]; then
    echo "This function takes one argument: a fastq filename" 1>&2
    usage
fi

if [ ! -f $1 ]; then
    echo "The fastq filename you provided does not exist" 1>&2
    usage
fi

## we print 1000 errors, should be enough
fastQValidator --noeof --file $1 --printableErrors 1000
