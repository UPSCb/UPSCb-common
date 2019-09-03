#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH -p core -n 1
#SBATCH -t 12:00:00

module load bioinfo-tools Reaper
set -ex

infile=$1
dir=$2

if [ ! -f $infile ]; then
    echo "invalid file"
    exit 1
fi

if [ ! -d $dir ]; then
    echo "invalid directory"
    exit 1
fi

name=$(basename $infile)

minion search-adapter -i $infile > $dir/${name/.fastq.gz/.adapter_mn.txt}
