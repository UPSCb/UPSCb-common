#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL

####
#	A runner part of the novel genen and long non coding RNA pipeline
#	python get_fasta_seq.py <gff_file> <fasta_file> <output>
####

## stop on error but be verbose
set -e
set -x

python $UPSCb/src/python/novel_genes/get_fasta_seq.py $1 $2 $3