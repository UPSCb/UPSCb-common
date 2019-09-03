#!/bin/bash -l

#SBATCH -p core -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e

usage() {
    echo >&2 "usage: $0 [OPTIONS] input.gff

Compress a gff, vcf, bed, sam or psltab file with bgzip
and index it with tabix.

OPTIONS:
-h          Show this message and exit.
-p FORMAT   Format of input. Valid arguments of FORMAT are
            gff, vcf, bed, sam and psltab. Default is gff."
}

if ! hash bgzip 2>/dev/null; then
    echo >&2 "ERROR: bgzip not found in path"
    exit 1
fi

if ! hash tabix 2>/dev/null; then
    echo >&2 "ERROR: tabix not found in path"
    exit 1
fi

OPTIND=1
format="gff"
while getopts "hp:" opt; do
    case "$opt" in
        h) usage; exit 1 ;;
        p) format=$OPTARG ;;
        ?) usage; exit 1 ;;
    esac
done

shift $((OPTIND - 1))

if [ $# -lt 1 ]; then
    echo >&2 "ERROR: input file missing"
    usage
    exit 1
fi

if [ ! -f $1 ]; then
    echo >&2 "ERROR: could not find input file: '$1'"
    usage
    exit 1
fi

bgzip $1
tabix -p $format $1.gz
