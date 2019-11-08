#!/bin/bash
#SBATCH -p core -n 8
#SBATCH -t 3:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL

set -ex

#module load bioinfo-tools kallisto samtools

## a usage function
usage(){
    echo >&2 \
"
    Usage: 
    Paired end: $0 [options] <fwd fastq file> <rev fastq file> <inx file> <transcript fasta file> <out dir>
    Single end: $0 [options] -s <fastq file> <inx file> <transcript fasta file> <out dir>
    Options:
      b: the number of boostrap; default to 100
      c: convert the BAM to CRAM (requires -p). Default FALSE
      f: force; i.e. overwrite results
      m: the memory multiplier, in GB; i.e. total memory = -m * -t (7 * 8 CPU = 56G by default )
      s: single end sequencing; the function only expects 3 arguments then
      t: the number of threads; default to 8
      p: generate pseudobam; default to FALSE
      u: unstranded data; default to Illumina rf-stranded
      r: Strand specific reads, first read reverse, corresponds to Illumina rf-stranded
      F: Strand specific reads, first read forward, corresponds to non Illunmina fr-stranded
      M: fragment length mean
      S: fragment length sd
" 
    exit 1
}

## VARS
BOOTSTRAP="-b 100"
CPU=8
CRAM=0
MEMMULT=8
SINGLE=0
FORCE=0
STRANDED="--rf-stranded"
THREADS="-t $CPU"
PSEUDOBAM=
FRAGMENT_LENGTH_MEAN=
FRAGMENT_LENGTH_SD=

# get the options
while getopts b:cfm:st:purFM:S: option
  do
    case "$option" in
    b) BOOTSTRAP="-b $OPTARG";;
    c) CRAM=1;;
    f) FORCE=1;;
    m) MEMMULT=$OPTARG;;
  	s) SINGLE=1;;
  	t) THREADS="-t $OPTARG";;
  	p) PSEUDOBAM="--pseudobam";;
  	u) STRANDED=;;
  	r) STRANDED="--rf-stranded";;
  	F) STRANDED="--fr-stranded";;
  	M) FRAGMENT_LENGTH_MEAN=$OPTARG;;
  	S) FRAGMENT_LENGTH_SD=$OPTARG;;
	  *) ## unknown flag
	    usage;;
  esac
done
shift `expr $OPTIND - 1`

## we get 3 files and one dir as input 
if [ $# -lt 4 ] || [ $# -gt 5 ] ; then
    echo "This function takes one (two) fastq file(s), one index file, one fasta file and one output dir as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing fastq file"
    usage
fi

if [ $SINGLE -eq 0 ]; then
  if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing fastq file"
    usage
  fi
  if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing index file"
    usage
  fi
  if [ ! -f $4 ]; then
    echo "The forth argument needs to be an existing fasta file"
    usage
  fi
  if [ ! -d $5 ]; then
    echo "The fifth argument needs to be an existing directory"
    usage
  fi
  
  fnam=$(basename ${1/_1.f*q.gz/})
  outdir=$5/$fnam
  if [ ! -d $outdir ]; then
    mkdir -p $outdir
  fi
  
  if [ ! -f $outdir/abundance.tsv ] || [ $FORCE -eq 1 ]; then
  
    kallisto quant -i $3 $BOOTSTRAP \
    -o $outdir $PSEUDOBAM $THREADS $STRANDED $1 $2 \
    > $outdir/$fnam.sam

    samtools view -bT $4 $outdir/$fnam.sam | \
    samtools sort -o $outdir/${fnam}_pseudo.bam -@ $CPU -m ${MEMMULT}G && \
    samtools index $outdir/${fnam}_pseudo.bam
  fi

else
  if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing index file"
    usage
  fi
  if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing fasta file"
    usage
  fi
  if [ ! -d $4 ]; then
    echo "The forth argument needs to be an existing directory"
    usage
  fi
  if [ -z $FRAGMENT_LENGTH_MEAN ] || [ -z $FRAGMENT_LENGTH_SD ]; then
    echo "if single-end reads you need to provide fragment length mean and sd using -M and -S"
    usage
  fi  
  
  fnam=$(basename ${1/.f*q.gz/})
  outdir=$4/$fnam
  if [ ! -d $outdir ]; then
    mkdir -p $outdir
  fi
  
  if [ ! -f $outdir/abundance.tsv ] || [ $FORCE -eq 1 ]; then
    kallisto quant -i $2 $BOOTSTRAP \
    -o $outdir $PSEUDOBAM $THREADS $STRANDED \
    --single $1 -l $FRAGMENT_LENGTH_MEAN -s $FRAGMENT_LENGTH_SD

    if [ ! -z $PSEUDOBAM ] && [ $CRAM -eq 1 ]; then  
      samtools view -CT $3 $outdir/pseudoalignments.bam | \
      samtools sort -o $outdir/${fnam}_pseudo.cram -@ $CPU -m ${MEMMULT}G && \
      samtools index $outdir/${fnam}_pseudo.cram
    
      if [ -f ${fnam}_pseudo.cram ]; then
        rm $outdir/pseudoalignments.bam
      fi
    fi
  fi
fi
