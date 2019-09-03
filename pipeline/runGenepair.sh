#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=ALL

set -e

print_help() {
    echo >&2 "usage: $0 [OPTIONS] expression genes samples

    Compute MI/CLR from a gene expression matrix

OPTIONS:

    -i      expression input is MI lower triangular matrix
            from previous calculation, compute CLR
    -b      number of bins for spline calculation (default: 7)
    -o      spline order (default: 3)
    -m      Only compute MI
    -c      Complete, i.e. no/few missing values

NOTES:

    The expression matrix must be space separated with missing
    values coded as 1e-30 and must not have header or row names.
    Each columns should correspond to a sample, and each row to
    a gene.

    The 'genes' and 'samples' files should be one-column text-files
    containing gene and sample names in the same order as the
    expression matrix.

    You must have the genepair executable in your path."
    exit 1
}

OPTIND=1
OPTARG="null"
mimethod=6
clr=true
bins=7
splineorder=3
ismi=false
while getopts "hmcb:o:i" opt; do
    case "$opt" in
        h) print_help ;;
        m) clr=false ;;
        c) mimethod=1 ;;
        b) bins=$OPTARG ;;
        o) splineorder=$OPTARG ;;
        i) ismi=true ;;
        ?) print_help ;;
    esac
done

shift $((OPTIND-1))

if [ $# -lt 3 ]; then
    echo >&2 "$0: error: too few arguments"
    print_help
fi

if ! hash genepair 2>/dev/null; then
    echo >&2 "$0: error: you don't have genepair in your path"
    print_help
fi

expression=$1
genes=$2
samples=$3
outdir=$(dirname "$expression")
fname=$(basename "${expression%.txt}")

[[ ! -f "$expression" ]] && echo >&2 "$0: error: expression file not found" && print_help
[[ ! -f "$genes" ]] && echo >&2 "$0: error: gene file not found" && print_help
[[ ! -f "$samples" ]] && echo >&2 "$0: error: sample file not found" && print_help

# If input is MI we can only compute CLR
if $ismi && ! $clr; then
    echo >&2 "$0: error: if input is MI, you must compute CLR"
    print_help
fi

# Test if bins and spline order are positive integers
re='^[0-9]+$'
if ! [[ $bins =~ $re ]]; then
    echo >&2 "$0: error: bins must be positive integer"
    exit 1
fi

if ! [[ $splineorder =~ $re ]]; then
    echo >&2 "$0: error: spline order must be positive integer"
    exit 1
fi

# Number of genes and samples
ngene=$(wc -l "$genes" | cut -f1 -d" ")
nsample=$(wc -l "$samples" | cut -f1 -d" ")

# MI calculation
if ! $ismi; then
    echo >&2 "command line:" genepair $mimethod "$expression" "${outdir}/${fname}_mi.txt" $ngene $nsample 0 $ngene $bins $splineorder
    genepair $mimethod "$expression" "${outdir}/${fname}_mi.txt" $ngene $nsample 0 $ngene $bins $splineorder
fi

if $ismi; then
    mifile=$expression
else
    mifile="${outdir}/${fname}_mi.txt"
fi

# CLR calculation
if $clr; then
    echo >&2 "command line:" genepair 2 "$mifile" $ngene "${mifile%.txt}_clr.txt" 1
    genepair 2 "$mifile" $ngene "${mifile%.txt}_clr.txt" 1
fi
