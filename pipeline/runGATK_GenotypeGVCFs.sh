#!/bin/bash
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL

# check -u
set -ex

# module load bioinfo-tools GATK

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

PLOIDY=2

USAGETXT=\
"
Usage: $0 [options] <ref.fa> <out.vcf> <combined gvcf>

Notes:This script is GATK v4 compatible and GATK V3 incompatible. 
This script was changed to use a CombinedGVCFs - check runGATK_CombineGVCFs.sh

Options: -p ploidy defaults to 1

"

## get the options
while getopts p: option
do
  case "$option" in
	    p) PLOIDY=$OPTARG;;
		\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

# check
isExec gatk

if [ "$#" -ne "3" ]; then
    usage
fi

if [ ! -f $1 ]; then
    abort "ERROR: could not find reference: '$1'"
fi

ref=$1
shift

out=$1
shift

in=$1
shift

# variants=()
# for gvcf in $@; do
#     if [ ! -f "$gvcf" ]; then
#         abort "ERROR: file not found: '$gvcf'"
#     fi
#     variants+=("--variant $gvcf")
# done

# checth GVCF options
gatk GenotypeGVCFs -R "$ref" -V $in -O $out --sample-ploidy $PLOIDY

