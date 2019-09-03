#!/bin/bash

#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL

usage() {
    echo >&2 "usage: $(basename $0) [options] vcf [vcf ...]

Print a markdown formatted table of SNP statistics for
one or more VCF files. This requires bcftools v0.2 or higher.

OPTIONS:
    -f      variant must have \"PASS\" or \".\" in its filter field

NOTE:
    If you are running this on UPPMAX, note that loading
    samtools 0.1.9 or higher will put an older version of
    bcftools in your path.

    If you are running this through SLURM, be sure to set
    the output to something meaningful since this is where
    the results will end up."
    exit 1
}

filter=
OPTIND=1
while getopts "f" opt; do
    case "$opt" in
        f) filter="-f PASS,.";;
        ?) usage;;
    esac
done

shift $((OPTIND - 1))

if ! hash bcftools 2>/dev/null; then
    echo >&2 "ERROR: bcftools not found"
    usage
fi

if ! bcftools --version 2>&1 | grep 'bcftools 0.2' >/dev/null 2>&1; then
    echo >&2 "ERROR: version 0.2+ of bcftools required"
    usage
fi

if [ $# -lt 1 ]; then
    echo >&2 "ERROR: at least one vcf file must be given"
    usage
fi

printf '| Sample | # SNPs | # MNPs | # indels | # others | # multiallelic sites | # multiallelic SNP sites |\n'
printf '|--------|--------|--------|----------|----------|----------------------|--------------------------|'
for f in $@; do
    printf "file:\t$f\n"
    bcftools stats $filter $f | grep "^SN" | tail -n +2 | cut -f3,4
done | \
awk 'BEGIN {FS = "\t"}; {if ($1 == "file:") {smpl=$2; sub(/^.*\//, "", smpl); sub(/.raw.*$/, "" smpl); printf("\n| %s |", smpl)} else {printf(" %s |", $2)}}' | \
sort -k1 -t'|'
