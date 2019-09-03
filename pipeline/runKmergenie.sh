#!/bin/bash -l
#SBATCH -c 8

module load bioinfo-tools kmergenie

usage ()
{
	echo "runKmergenie.sh <READ1> [READ2...] <OUT_DIR>" >&2
	echo "" >&2
	exit 1
}

NARGS=$#
CALL=$@

if [ $NARGS -lt 2 ]; then
	echo -e "\e[91m[ERROR] The minimum number of arguments is 2\e[39m" >&2
	usage
fi

i=0
let NARGS-=1
while [ $i -lt $NARGS ]; do
	readFilesIn[$i]=$1
	shift
	let i+=1
done

OUT_DIR=$@

if [ ! -d $OUT_DIR ];then
	echo -e "\e[33m[WARN] $OUT_DIR not a directory... creating\e[39m" >&2
	mkdir -p $OUT_DIR
fi

for f in ${readFilesIn[@]}; do
	if [ ! -f $f ];then
		echo -e "\e[91m[ERROR] \"$f\" is not a valid file\e[39m" >&2
		usage
	fi
done

TMF=$(tempfile)
echo $TMF
for f in ${readFilesIn[@]}; do
	echo -e $(readlink -f $f) >> $TMF
done
BN=$(basename  ${readFilesIn[0]})

echo -e "[INFO] Shell:	$0 $CALL" >&2
echo -e "[INFO] Exec: kmergenie --diploid -t 8 $TMF -o ${OUT_DIR}/${BN}_hist" >&2
echo -e "[INFO] PWD:  	$PWD" >&2
kmergenie --diploid -t 8 $TMF -o ${OUT_DIR}/${BN}_hist
rm $TMF
