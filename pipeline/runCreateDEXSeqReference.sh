#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.sundell@umu.se

module load bioinfo-tools
module load python/2.6.6

## abort on error
set -e

## usage

echo "runDEXSeq_count.sh input_gff output_gff"

python /mnt/picea/home/ishutava/Git/UPSCb/src/python/dexseq_prepare_annotation.py /mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/gff/TAIR10_GFF3_genes_transposons.gtf /mnt/picea/projects/docker/upsc2017/jBrowse/srobert/dr4-resistant-mutant/TAIR10_GFF3_MY.gff

