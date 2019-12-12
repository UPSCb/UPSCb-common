#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mem 6G

set -e

#module load bioinfo-tools GATK
#module load java

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# Defaults
JavaMem=6G
USAGETXT=\
"
Usage: $0 <BAM file> <reference fasta> <target interval> <output directory>

Note: This script is not GATK v4 compatible. Load a GATK V3 module. More at https://software.broadinstitute.org/gatk/blog?id=7847
"

# GATK 3
if [ -z $GATK_HOME ]; then
  usage
fi
GATK=$GATK_HOME/GenomeAnalysisTK.jar

if [ -d "$SNIC_TMP" ]; then
    tmp=$SNIC_TMP
else
    tmp=/mnt/picea/tmp
fi

if [ $# -lt 4 ]; then
  usage
fi

if [ ! -f $1 ]; then
    abort "Could not find BAM file '$1'"
fi
inbam=$1

if [ ! -f $2 ]; then
    abort "Could not find reference '$2'"
fi
ref=$2

if [ ! -f $3 ]; then
    abort "Could not find interval file '$3'"
fi
interval=$3


if [ ! -d $4 ]; then
    abort "Could not find directory '$4'"
fi
outdir=$4

# drop all four args 
shift
shift
shift
shift

inname=${inbam##*/}
outname=${inname%.bam}

# Perform local realignment
java -Xmx${JavaMem} -jar -Djava.io.tmpdir=$tmp $GATK \
    -T IndelRealigner \
    -I $inbam \
    -R $ref \
    --targetIntervals $interval \
    -o $outdir/${outname}_realigned.bam $@

