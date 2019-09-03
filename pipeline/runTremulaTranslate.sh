#!/bin/bash -l

set -eu

#SBATCH -t 40:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL

if [ $# -lt 3 ]; then
    echo >&2 "usage: $0 <bam> <transtable> <out.bam>"
    exit 1
fi

inbam=$1
transtable=$2
outbam=$3

[[ ! -f $inbam ]] && echo "Could not find BAM file" && exit 1
[[ ! -f $transtable ]] && echo "Could not find translation table" && exit 1
[[ ! -d $(dirname $outbam) ]] && echo "Could not find directory for output" && exit 1
[[ -z $UPSCb ]] && echo "The UPSCb environment variable needs to be set" && exit 1

module load bioinfo-tools samtools

perl $UPSCb/src/perl/tremula_scaffold_bamconvert.pl <(samtools view -h $inbam) $transtable | \
    samtools view -hSb - > $outbam
