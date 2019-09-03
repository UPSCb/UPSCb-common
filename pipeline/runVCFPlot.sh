#!/bin/bash

#SBATCH -p core -n 1
#SBATCH -t 2:00:00
#SBATCH --mail-type ALL

module load R/3.1.0

usage() {
    if [[ ! -z $1 ]]; then
        echo >&2 $1
    fi
    echo >&2 "usage: $0 input.vcf[.gz] outdir title"
    exit 1
}

[[ $# -lt 3 ]] && usage "ERROR: too few arguments"
[[ ! -f $1 ]] && usage "ERROR: file not found: $1"
[[ ! -d $2 ]] && usage "ERROR: output directory not found: $2"

Rscript ../../../src/R/plotVCFQual.R "$1" "$2" "$3"
