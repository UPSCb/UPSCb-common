#!/bin/bash -l

#SBATCH -t 1-00:00:00
#SBATCH -p core
#SBATCH -n 6
#SBATCH --mem 36G

set -ex

 module load java
 module load bioinfo-tools
 module load GATK
 module load R

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage: $0 <BAM file> <fasta ref> <output directory> <dbsnp>

Note: This script is  GATK v4 compatible and V3 incompatible. More at https://software.broadinstitute.org/gatk/blog?id=7847
"

# check GATK V4
isExec gatk

# Tests
if [ ! -f $1 ]; then
    abort "Could not find BAM file '$1'"
fi
inbam=$1

if [ ! -f $2 ]; then
    abort "Could not find FASTA file '$2'"
fi
ref=$2

if [ ! -d $3 ]; then
    abort "Could not find directory '$3'"
fi
outdir=$3

if [ ! -f $4 ]; then
    abort "Could not find dbSNP file '$4'"
fi
db=$4

name_out=`basename "${inbam/.bam/.table}"`
outname=`basename "${inbam/.bam/_recalibrated.bam}"`
post_out=`basename "${inbam/.bam/_recalibrated.table}"`
plots_out=`basename "${inbam/.bam/.pdf}"`

# Run BaseRecalibrator
## Analyze patterns of covariation in the sequence dataset
gatk BaseRecalibrator -R $ref --known-sites $db -I $inbam -O $outdir/$name_out

## Apply the recalibration to your sequence data
gatk ApplyBQSR -R $ref -I $inbam --bqsr-recal-file $outdir/$name_out -O $outdir/$outname

## Do a second pass to analyze covariation remaining after recalibration
gatk BaseRecalibrator -R $ref --known-sites $db -I $outdir/$outname -O $outdir/$post_out

## Generate before/after plots
gatk AnalyzeCovariates -before $outdir/$name_out -after $outdir/$post_out -plots $outdir/$plots_out
