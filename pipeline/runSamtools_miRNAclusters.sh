#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH -p core
#SBATCH -n 1

module load bioinfo-tools
module load samtools/1.3.1

file=$1
bed=$2
outdir=$3
name=$4

# extract alignments in miRNA loci regions
# extract columns with sequence names and genomic location
samtools view -L $bed $file | cut -f 1,3,4 > $outdir/$name.miRNA.txt
