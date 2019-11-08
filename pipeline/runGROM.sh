#!/bin/bash
#SBATCH -p core
#SBACTH -n 16

set -ex

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage: $0 <in bam> <fasta reference> <out dir>

Note: the GROM, htslib and vcftools modules need to be loaded
"

# Checks
isEnvVarSet("UPSCb")

isExec("GROM")

isExec("vcf-sort")

isExec("bgzip")

if [ $# != 3 ]; then
  abort "This script expects 2 arguments"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing bam file"
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be an existing fasta file"
fi

if [ ! -d $3 ]; then
  abort "The second argument needs to be an existing directory"
fi

# Setup
GROM=$(which GROM)
cd $out
ln -s $GROM.

# Run (P is half of the cores used)
outfile=$3/$(basename ${1/.bam/.vcf})
$GROM -M -P 8 -i $1 -r $2 -o $outfile

# Sort, compress and index
sorted=${outfile/.vcf/_sorted.vcf}
cat $outfile | vcf-sort -c > $sorted

bgzip -f $sorted

tabix ${sorted}.gz -p vcf
 
# Cleanup
find $out -type l -name "GROM" -delete
rm $outfile
