#!/bin/bash -l
#SBATCH -n 1


set -eu

# You need a singularity container for blast, a genome fasta file,
# a fasta file with the sequence of the used gRNAs and the name of an output file
container=$(realpath SINGULARITY_CONTAINER)
genome=$(realpath GENOME_FASTA_FILE)
gRNA=$(realapth gRNA_FASTA_FILE)
output="NAME_OUTPUT_FILE"

apptainer exec -B /mnt $container makeblastdb -in $genome -dbtype nucl

apptainer exec -B /mnt $container blastn -query  $gRNA \
-db $genome -task blastn-short -outfmt 7 -out $output
