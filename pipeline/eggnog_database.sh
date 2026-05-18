#!/bin/bash
#SBATCH -A account
#SBATCH -t 48:00:00
#SBATCH -o download_eggnog.out
#SBATCH -e download_eggnog.err



ml GCC/13.2.0 OpenMPI/4.1.6

ml eggnog-mapper/2.1.12

download_eggnog_data.py -y -P -F -M -H -d DATABASE  --data_dir data/reference/eggnog_database/
