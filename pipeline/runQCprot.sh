#!/bin/bash -l
#SBATCH -p all
#SBATCH -n 8 --mem=150G -t 2-00:00:00
## no mail at the time
###SBATCH --mail-type=ALL

## stop on error
set -e

## be verbose and extend the commands
set -x

## load the modules
module load bioinfo-tools
module load blast

blastx -query $1/Trinity.format.fasta \
         -db $2 -out $3/blastxprot.outfmt6 \
         -evalue 1e-20 -num_threads 8 -max_target_seqs 1 -outfmt 6

