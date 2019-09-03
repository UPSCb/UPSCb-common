#!/bin/bash -l
#SBATCH -c 8
module load bioinfo-tools ABySS

usage () 
{
    echo "runAbyssBloom.sh <k> <out file> <read1> <read2> ..."
    echo
    exit 1
}

K=$1
OUT=$2

if [ ! $K =~ '^[0-9]+$' ]; then
    echo "First argument has to be a number"
    usage
fi

shift 2

for f in $@; do
    if [ ! -f $f ];then
	echo "$f is not a valid file"
	usage
    fi
done

if [ ! -d $(dirname $OUT) ];then
    echo "Creating out directory"
    mkdir -p $(dirname $OUT)
fi

abyss-bloom build -k $K -t 8 $OUT $@


