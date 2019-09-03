#!/bin/bash -l
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH -c 1

set -ex

module load perl bioinfo-tools sspace

usage()
{
	echo "Usage:$0 <contigs file> <out dir> <gzipped pacbio reads>"
	exit 1
}

CONTIGS=$1
OUT=$2
shift 2
PB=( "$@" )

for f in ${PB[@]}; do
	if [ ! -f $f ]; then
	echo "$f not a valid file"
	exit 1
	fi
done

if [ ! -f $CONTIGS ]; then
	echo "Contigs file $CONTIGS does not exist"
	exit 1
fi

if [ ! -d $OUT ]; then
	echo "Out dir $OUT does note exist"
	exit 1
fi

F1=${PB[0]}
EXT="${F1##*.}"

if [[ $EXT == "gz" ]]; then
	EXT=${F1/.gz/}
	EXT=${EXT##*.}
	GZIP=1
else
	GZIP=0
fi

touch $OUT/tmp.ct.$EXT
if [ $GZIP -eq 1 ];then
	zcat ${PB[@]} > $OUT/tmp.ct.$EXT
else
	cat ${PB[@]} > $OUT/tmp.ct.$EXT
fi

perl $(which SSPACE-LongRead.pl) -c $CONTIGS -p $OUT/tmp.ct.$EXT -b $OUT

rm $OUT/tmp.ct.$EXT
