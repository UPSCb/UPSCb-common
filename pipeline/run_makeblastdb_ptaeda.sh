#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools
module load blast/2.2.29+

# run the command
makeblastdb -in /mnt/picea/storage/reference/Pinus-taeda/v1.01/fasta/ptaeda.v1.01-genome-collapsed-for-STAR.fa -dbtype nucl -out /mnt/picea/projects/spruce/pipeline/psari_data/pinus_taeda/db_blastplus/ptaeda.v1.01-genome.fa


