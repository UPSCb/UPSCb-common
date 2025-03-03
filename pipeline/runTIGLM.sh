#!/bin/bash -l
#SBATCH -p main
#SBATCH -c 8
#SBATCH --mail-type=ALL

## stop on error
set -ex

## we get one dir and one file as input
usage(){
  echo >&2 \
  " Usage: $0 [options] <expression matrix file> <gene names file> <out dir>
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
if [ "$#" -lt 3 ]; then
  echo "This function requires 3 arguments"
  usage;
fi

if [ ! -f $1 ]; then
  echo "The first argument needs to be an existing file"
  usage;
fi

if [ ! -f $2 ]; then
  echo "The second argument needs to be an existing file"
  usage;
fi

if [ ! -d $3 ]; then
  echo "The third argument needs to be an existing directory"
  usage;
fi

# set up
cd $3
mkdir tmp

# run
/mnt/picea/home/bastian/Git/geneNetworkR/src/tiglm/bin/tiglm -i $1 -g $2 -p $CPU

# clean up
mv MA* tmp

# and concatenate
cat tmp/MA* > TIGLM-edgelist.tsv

# finally clean up
rm -rf tmp/MA*
rmdir tmp


