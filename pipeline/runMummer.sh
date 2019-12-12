#!/bin/bash -l
#SBATCH -n 1 -p core
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=16G


## be verbose and stop on error
set -ex

## modes
modes=(mummer mummer1 mummer3 mummer3D nucmer)

## usage
export USAGETXT=\
"
Usage: $0 <mode> <query fasta> <reference fasta> <out dir>

<mode> is one of
${modes[@]}

Note:
   You need to set the UPSCb env. variable to your UPSCb git checkout directory
"

## common functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## arguments
if [ $# != 4 ]; then
    abort "This function takes 3 arguments: the mode, the query fasta and the reference fasta files and the output directory"
fi

MODE=$1
shift
if [ $(containsElement $MODE "${modes[@]}") -eq 1 ]; then
  abort "Unknown mode"
fi

QRY=$1
shift
if [ ! -f $QRY ]; then
  abort "The second argument must be a valid fasta file"
fi

REF=$1
shift
if [ ! -f $REF ]; then
  abort "The third argument must be a valid fasta file"
fi

OUT=$1
shift
if [ ! -d $OUT ]; then
    abort "The third argument must an existing directory"
fi

# check for compressed query
QRYext="${QRY##*.}"
if [ "$QRYext" == "gz" ]; then
  ## create temporary file
  QRYtmp=`tempfile`

  ## decompress
  gunzip -c $QRY > $QRYtmp
else
  QRYtmp=$QRY
fi

# check for compressed ref
REFext="${REF##*.}"
if [ "$REFext" == "gz" ]; then
  ## create temporary file
  REFtmp=`tempfile`

  ## decompress
  gunzip -c $REF > $REFtmp
else
  REFtmp=$REF
fi

cd $OUT

## command
case "$MODE" in
  "mummer")
    mummer -mum -b -c $REFtmp $QRYtmp > ref_qry.mums
    mummerplot --postscript --prefix=ref_qry ref_qry.mums
    gnuplot ref_qry.gp
  ;;
  "mummer1")
    run-mummer1 $REFtmp $QRYtmp ref_qry
    run-mummer1 $REFtmp $QRYtmp ref_qry_rev -r
  ;;
  "mummer3")
    run-mummer3 $REFtmp $QRYtmp ref_qry
  ;;
  "mummer3D")
    run-mummer3D $REFtmp $QRYtmp ref_qry
  ;;
  "nucmer")
    nucmer --maxgap=500 --mincluster=100 --prefix=$OUT/ref_qry $REFtmp $QRYtmp
    show-coords -r $OUT/ref_qry.delta > $OUT/ref_qry.coords
    #show-aligns ref_qry.delta refname qryname > ref_qry.aligns
    #delta-filter -q -r ref_qry.delta > ref_qry.filter
    #mummerplot ref_qry.filter -R ref.fasta -Q qry.fasta
  ;;
esac
