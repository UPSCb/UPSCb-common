#!/bin/bash -l
#SBATCH -p core
#SBATCH -c 16
#SBATCH --mail-type=ALL

## stop on error
set -ex

## we get input and output file as input
usage(){
  echo >&2 \
  " Usage: $0 [options] <expression matrix file> <out file>
    Options: -c set the number of CPUs
  "
  exit 1
}

## Set defaults
CPU=16

## get the options
while getopts c: option
  do
    case "$option" in
	    c) CPU=$OPTARG;;
		  \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## Process the args
if [ "$#" -lt 2 ]; then
  echo "This function requires 3 arguments"
  usage;
fi

if [ ! -f $1 ]; then
  echo "The first argument needs to be an existing file"    
  usage;
fi

## The R script now checks output permission so no
## need to do it here

## check the environment
if [ -z $UPSCb ]; then
  echo "You need to define the UPSCb env var to your local UPSCb git checkout dir"
  usage
fi

## load the modules
module load R

## and run
Rscript --vanilla $UPSCb/src/R/GENIE3_simple.R $1 $2 $CPU

## cleanup
cat tmp/file* > GENIE3-edgelist.tsv
rm -rf tmp

