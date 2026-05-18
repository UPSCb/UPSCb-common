#!/bin/bash
#SBATCH -A account
#SBATCH -t 168:00:00
#SBATCH --output eggnog_log.out
#SBATCH --error eggnog_log.err
#SBATCH -c 20




ml GCC/13.2.0 OpenMPI/4.1.6

ml eggnog-mapper/2.1.12

emapper.py --temp_dir temp --data_dir eggnog_database/ --genepred search  \
 -i  proteome.fa \
--cpu 20 -o eggnog_annotation --output_dir outdir -m diamond --itype proteins --override

