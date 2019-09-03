#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.sundell@umu.se

module load bioinfo-tools
module load HTSeq/0.6.1

## abort on error
set -e

## usage

echo "runDEXSeq_count.sh input_file input_gff"

#python ~/script/python/dexseq_count.py $1 $2 -p yes -f bam -s no

name1=${1##*/}
name=${name1%.bam}


#echo $3/$name.txt
python ~/Git/UPSCb/src/python/dexseq_count.py $2 $1 $3/$name".txt" -p yes -f bam -s no
