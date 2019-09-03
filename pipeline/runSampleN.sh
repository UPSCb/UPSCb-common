#!/bin/bash

usage(){
  echo "$0 <in> <out> <file prefix> <subset size in M reads>"
  exit 1
} 

# sanity check
if [ $# -ne 4 ]; then
  echo "This script expects 4 arguments"
  usage
fi

# the vars
in=$1
out=$2
nam=$3
subset=$4

# sanity check
if [ ! -d $in ]; then
  echo "The first arg needs to be the input dir"
  usage
fi

if [ ! -d $out ]; then
  echo "The second arg needs to be the output dir"
  usage
fi

if [ ! -f $in/${nam}_1.fq.gz ]; then
  echo "The third arg needs to be the prefix (without _1.fq.gz) of the input files"
  usage
fi

# run
sampleN -l $out/${nam}_${subset}.log.txt.gz -n `expr $subset "*" 1000000` -o $out/${nam}_${subset}_million $in/${nam}_1.fq.gz $in/${nam}_2.fq.gz 

