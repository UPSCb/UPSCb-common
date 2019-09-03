#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH -p all
#SBATCH -n 8
#SBATCH -t 2-00:00:00

# module load bioinfo-tools Reaper 
set -ex

meme $1.$SLURM_ARRAY_TASK_ID -dna -oc $2$SLURM_ARRAY_TASK_ID -mod anr -evt 0.05 -maxsize 3500000 -maxw 30 -nmotifs 100 -bfile $3 -p 8
