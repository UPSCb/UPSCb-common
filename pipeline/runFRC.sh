#!/bin/bash -l

## report error
set -e

## be verbose 
set -x

## usage
usage(){
echo >&2 \
"
     Usage: runFRC.sh <alignment bam> <out dir>
     Note: At the moment it only accept --pe-sam
"
exit 1
}

## check params
if [ $# != 2 ]; then
    echo "This script needs two parameters"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first parameter should be a bam file"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument has to be an existing output directory"
    usage
fi

# load modules
module load bioinfo-tools FRC

## run
cd $2
FRC --pe-sam $1 \
--pe-max-insert 425 --genome-size 11350000 --CEstats-PE-min -5 \
--CEstats-PE-max 5.5 --output i425g1135Cmin5Cmx55



