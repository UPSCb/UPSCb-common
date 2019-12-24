#!/bin/bash -l

#SBATCH -t 1-00:00:00
#SBATCH -p core
#SBATCH -n 6
#SBATCH --mem 36G

set -e

#module load java
#module load bioinfo-tools
#module load GATK

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: $0 <BAM file> <fasta ref> <output directory>

Note: This script is not GATK v4 compatible. Load a GATK V3 module. More at https://software.broadinstitute.org/gatk/blog?id=7847
"

# default
Threads=6
JavaThreadMem=6G

# GATK 3
if [ -z $GATK_HOME ]; then
  usage
fi
GATK=$GATK_HOME/GenomeAnalysisTK.jar

if [ $# -lt 3 ]; then
  usage
fi

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

# drop the three args
shift
shift
shift

name_out=`basename "${inbam/.bam/.intervals}"`

# Run
java -Xmx${JavaThreadMem} -jar $GATK -nt $Threads -I $inbam -R $ref -T RealignerTargetCreator -o $outdir/$name_out $@

