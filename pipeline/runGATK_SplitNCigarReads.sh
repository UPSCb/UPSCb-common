#!/bin/bash -l

#SBATCH -t 1-00:00:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem=32G

set -e

module load bioinfo-tools GATK
module load java

if [ -z $GATK_HOME ]; then
    echo >&2 "Could not find GATK"
    exit 1
fi

GATK=$GATK_HOME/GenomeAnalysisTK.jar
if [ -d "$SNIC_TMP" ]; then
    tmp=$SNIC_TMP
else
    tmp=/mnt/picea/tmp
fi

if [ $# -ne 3 ]; then
    echo "Usage: $0 <in.bam> <ref.fasta> <output directory>" 1>&2
    exit 1
fi

if [ ! -f $1 ]; then
    echo "Could not find BAM file '$1'" 1>&2
    exit 1
fi
inbam=$1

if [ ! -f $2 ]; then
    echo "Could not find reference '$2'" 1>&2
    exit 1
fi
ref=$2

if [ ! -d $3 ]; then
    echo "No such output directory" 1>&2
    exit 1
fi
outdir=$3

namein=${inbam/.bam/}
bname=`basename $inbam`
sname=${bname/_[st]*[st]*_STAR*.bam/}

# Run splitNCigarReads
java -jar -Xmx5G -Djava.io.tmpdir=$tmp $GATK -T SplitNCigarReads \
    -R "$ref" \
    -I "$inbam" \
    -o "$outdir/${bname/.bam/}_split.bam" \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS
