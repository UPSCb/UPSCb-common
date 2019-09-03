#!/bin/bash

# usage 
USAGETXT=\
"
runSampleN_SE <out> <file> <subset size in M reads>

Note:
  Adapted to SE by mquevedo
"
 

# sanity check
if [ $# -ne 3 ]; then
  echo "This script expects 3 arguments"
  usage
fi

# the vars

out=$1
nam=$2
fnam=$(basename ${nam/_trim.fq.gz/}) 
subset=$3

echo
# sanity check
if [ ! -d $out ]; then
  echo "The first arg needs to be the output dir"
  usage
fi

if [ ! -f $nam ]; then
  echo "The second arg needs to be the file"
  usage
fi

# run
sampleN -n `expr $subset "*" 1000000` -o $out/$fnam"_"$subset"M" $nam 

