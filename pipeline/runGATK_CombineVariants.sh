#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL

if [ -z $SLURM_SUBMIT_DIR ]; then
    ## Not using SLURM
    if [ $# -lt 5 ]; then
        echo >&2 "Usage: $0 <gatk> <tmp> <ref.fasta> <out.vcf> <in.vcf> <in.vcf> [<in.vcf> ...]"
        exit 1
    fi
    GATK=$1
    tmp=$2
    shift 2
else
    if [ $# -lt 3 ]; then
        echo "Usage: $0 <ref.fasta> <out.vcf> <in.vcf> <in.vcf> [<in.vcf> ...]"
        exit 1
    fi
    GATK="/sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar"
    tmp=$SNIC_TMP
fi

if [ ! -f $GATK ]; then
    echo >&2 "ERROR: could not find GATK"
    exit 1
fi

if [ ! -d $tmp ]; then
    echo >&2 "ERROR: tmp is not a directory"
    exit 1
fi

if [ ! -f $1 ]; then
    echo >&2 "ERROR: could not find reference: '$1'"
    exit 1
fi
ref=$1

if [ ! -d `dirname $2` ]; then
    echo >&2 "ERROR: no such directory: '$2'"
    exit 1
fi
out=$2

shift 2

variants=()
for f in $@; do
    variants=("${variants[@]}" "-V $f")
done

VAR=(${variants[@]})

java -Xmx4g -jar $GATK \
    -T CombineVariants \
    -R $ref \
    --filteredAreUncalled \
    "${VAR[@]}" \
    -o $out
