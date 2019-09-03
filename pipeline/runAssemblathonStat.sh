#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## load the modules
module load perl

## we get one dir and one file as input
usage(){
    echo >&2 \
    " Usage: $0 [options] <genome fasta file>
      
      Options:
        -n the maximum number of consecutive N characters allowed 
           before scaffolds are split into contigs (default 25)
        -g the estimated genome size
        
      Notes: The result will be written in the same dir as the input fasta file  
    "
    exit 1
}

OPTIONS="-csv"
# getting the options
while getopts n:g: option
do
        case "$option" in
      g) OPTIONS="$OPTIONS -genome_size $OPTARG";;
	    n) OPTIONS="$OPTIONS -n $OPTARG";;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

if [ "$#" -lt 1 ]; then
  echo "This function requires 2 arguments"
  usage;
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing fasta file"    
    usage;
fi

## get the subtracted results
perl -I $UPSCb/src/perl $UPSCb/src/perl/assemblathon_stats.pl $OPTIONS $1
