#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH -p core -n 1
#SBATCH -t 2-00:00:00

# module load bioinfo-tools Reaper 
set -ex

infile=$1
dir=$2
shift 2

if [ ! -f $infile ]; then
    echo "invalid file"
    exit 1
fi

if [ ! -d $dir ]; then
    echo "invalid directory"
    exit 1
fi

name=$(basename $infile)

reaper -geom no-bc -i $infile -3pa TGGAATTCTCGGG -basename $dir/${name/.fastq.gz/} -nnn-check 1/1 -3pa ""
#reaper -i $infile -3pa TGGAATTCTCGGGTGCCAAGG -geom no-bc -basename $dir/${name/.fastq.gz/}  -nnn-check 1/1 -3pa "" -tabu TGGAATTCTCGGG $@
