#!/bin/bash -l 
#SBATCH -p nextflow
#SBATCH -t 72:00:00
#SBATCH -A SLURM_project_Code
#SBATCH -o log_rnaseq.out
#SBATCH -e log_rnaseq.err

set -eu -o pipefail

nextflow run nf-core/rnaseq -r 3.19.0 \
-profile singularity,upscb -c "nextflow/upscb.config" \
-params-file "nextflow/nf-params.json"  \
-with-trace -with-report "report_rnaseq.html" \
-resume
