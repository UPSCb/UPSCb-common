#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL

if [ -z $SLURM_SUBMIT_DIR ]; then
    if [ $# -lt 4 ]; then
        echo "Usage $0 </path/to/gatk> <tmp dir> <in.vcf> <ref.fasta> <output directory> [vf_args ...]" 1>&2
        exit 1
    fi

    GATK=$1
    tmp=$2
    shift 2
else
    if [ $# -lt 3 ]; then
        echo "Usage: $0 <in.vcf> <ref.fasta> <output directory> [vf_args ...]" 1>&2
        exit 1
    fi
    GATK=/sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar
    tmp=$SNIC_TMP
fi

if [ ! -d $tmp ]; then
    echo "tmp is not a directory" 1>&2
    exit 1
fi

if [ ! -f $1 ]; then
    echo "Could not find VCF file '$1'" 1>&2
    exit 1
fi
invcf=$1

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

shift 3

#vf_args="$@"

bname=`basename $invcf`
sname="${bname/.vcf*/}"
outfile="$outdir/${sname}_filtered.vcf"

# Run VariantFiltration
java -jar -Xmx2G -Djava.io.tmpdir=$tmp $GATK -T VariantFiltration \
    -V "$invcf" \
    -R "$ref" \
    -o "$outfile" \
    "$@" #${vf_args[@]}

