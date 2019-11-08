#!/bin/bash -l
#SBATCH -t UNLIMITED
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem 16G
#SBATCH --mail-type=ALL

set -e

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# default
JavaThreadMem=16G

# usage
export USAGETXT=\
"
Usage: $0 <VCF file> <interval file> <fasta ref> <output directory> [GATK additional options]

Note: several interval files can be provided, comma separated
"

if [ $# -lt 4 ]; then
  abort "This script expects 4 arguments"
fi

if [ ! -f $1 ]; then
    abort "Could not find VCF file '$1'"
fi
vcf=$1

# Check the second argument
echo $2 | xargs -d, -I {} bash -c 'if [ ! -f $0 ]; then abort "INTERVAL file $0 does not exist"; fi' {}
ivl=$2

if [ ! -f $3 ]; then
    abort "Could not find FASTA file '$3'"
fi
ref=$3

if [ ! -d $4 ]; then
    abort "Could not find directory '$4'"
fi
out=$4

# drop the four args
shift
shift
shift
shift

# run once for every interval file
for i in `echo $ivl | tr ',' ' '`; do

  nam=$out/$(basename $i)"-"$(basename "${vcf/.vcf/.fasta}")

  # Run
  java -Xmx${JavaThreadMem} -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T FastaAlternateReferenceMaker -R $ref -L $i -V $vcf  -o $nam $@
done
