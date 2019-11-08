#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL

set -eux

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: $0 <ref.fa> <out.vcf> <gvcf> [<gvcf> ...]

Notes:This script is GATK v4 compatible and GATK V3 incompatible.
"

# check
isExec gatk

if [ $# -lt 3 ]; then
    usage
fi

if [ ! -f $1 ]; then
    abort "ERROR: could not find reference: '$1'"
fi

ref=$1
shift

out=$1
shift

variants=()
for gvcf in $@; do
    if [ ! -f "$gvcf" ]; then
        abort "ERROR: file not found: '$gvcf'"
    fi
    variants+=("-V $gvcf")
done

# checth GVCF options
gatk CombineGVCFs -R "$ref" ${variants[@]} -O $out

