#!/bin/bash
#SBATCH --mail-type=all
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

set -ex

# module load bioinfo-tools java

usage(){
echo >&2 \
"
	Usage: $0 file config output
"
	exit 1
}

if [ $# -ne 3 ]; then
	echo "This file expects 3 arguments"
	usage
fi 

if [ ! -f $1 ]; then
	echo "The first argument should be an existing file"
	usage
fi
f=$1

if [ ! -f $cfg ]; then
	echo "The second argument should be an existing file"
	usage
fi
cfg=$2

if [ ! -d $3 ]; then
	echo "The third argument should be an existing dir"
	usage
fi
out=$3

name=$(basename ${f//.lane.clean.fa/})
java -Xmx8G -jar /mnt/picea/Modules/apps/bioinfo/srna-workbench/3.2/Workbench.jar -tool filter -srna_file $f -out_file $out/${name}_filtered.fa -params $cfg
