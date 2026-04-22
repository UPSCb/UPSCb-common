#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -t 168:00:00
#SBATCH --mem 2GB
#SBATCH -e nextflow_log.err
#SBATCH -o nextflow_log.out

set -u -e -o pipefail

ml Nextflow/25.10.0

nextflow run nf-core/rnaseq -r 3.24.0 \
  -profile hpc2n,singularity,gpu -params-file nextflow/rnaseq_STAR_GPU_potra_v2.json  -c nextflow/gpu.config -work-dir data/nobackup/workdir