#! /bin/env bash

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# DEFAULTS
CPU=16
trim="SLIDINGWINDOW:5:20 MINLEN:100"

# load the modules
# module load bioinfo-tools deML fastQC trimmomatic

# usage function
usage(){
echo >&2 \
"
  Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <I1.fastq.gz> <I2.fastq.gz> <inx file> <adapter file> <out dir>

"
	exit 1
}

## TODO we need a 16S/ITS switch and a link to references

# check the arguments
if [ $# -ne 7 ]; then
	usage
fi

if [ ! -f $1 ]; then
	echo "The forward fastq file: $1 does not exist"
	usage
fi

if [ ! -f $2 ]; then
  echo "The reverse fastq file: $2 does not exist"
	usage
fi

if [ ! -f $3 ]; then
  echo "The forward index file: $3 does not exist"
	usage
fi

if [ ! -f $4 ]; then
  echo "The reverse index file: $4 does not exist"
	usage
fi

if [ ! -f $5 ]; then
  echo "The mapping file: $5 does not exist"
	usage
fi

if [ ! -f $6 ]; then
  echo "The adapter file: $6 does not exist"
  usage
fi

if [ ! -d $7 ]; then
  echo "The output directory: $7 does not exist"
  usage
fi

# create out dirs
mkdir -p $7/FastQC/raw
mkdir -p $7/deML
mkdir -p $7/FastQC/deML
mkdir -p $7/trimmomatic
mkdir -p $7/FastQC/trimmomatic
mkdir -p $7/flash
mkdir -p $7/FastQC/flash
mkdir -p $7/OTU

# FastQC
fastqc -t 2 -o $7/FastQC/raw $1 $2

# rev-comp the I1 barcode
Rscript createMapFileForDeML.R -f $5 -o $7/deML

# deML
for f in `find $7/deML -type f -name "*_deML.txt"`; do
 deml/src/deML -i $f -f $1 -r $2 -if1 $3 -if2 $4 -o $7/deML/demultiplexed
 rm -f $7/deML/demultiplexed_unknown*
done

# fastqc
fastqc -t $CPU -o $7/FastQC/deML $7/deML/demultiplexed*_r[1,2]*.fq.gz

# trimmomatic
clip=ILLUMINACLIP:${6}:2:30:10
for f in `find $7/deML -type f -name "*_r[1,2].fq.gz"`;
do
  echo "${f//_r[1,2].fq.gz/}" ; done | sort | uniq | while read file;
  do
    pattern=`basename $file`

    java -jar Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $CPU -phred33 ${file}_r1.fq.gz \
    ${file}_r2.fq.gz $7/trimmomatic/${pattern}_trimmomatic_1.fq \
    $7/trimmomatic/${pattern}_unpaired_1.fq $7/trimmomatic/${pattern}_trimmomatic_2.fq \
    $7/trimmomatic/${pattern}_unpaired_2.fq $clip $trim

    find  $7/trimmomatic -type f -name "${pattern}_*.fq" | xargs -I {} -P $CPU gzip -f {}
  done

# fastqc
fastqc -t $CPU -o $7/FastQC/trimmomatic $7/trimmomatic/*_trimmomatic_[1,2].fq.gz

# flash
for f in $(find $7/trimmomatic -type f -name "*trimmomatic_[1,2].fq.gz"); 
  do echo ${f//_[1,2].fq.gz/}; done | sort | uniq | while read file; 
  do 
    pattern=`basename $file`_flash 
    flash -t $CPU -m 30 -M 251 -o $pattern -d $7/flash -O \
    ${file}_1.fq.gz ${file}_2.fq.gz > $7/flash/$pattern.log
    zcat $7/flash/$pattern.notCombined_1.fastq.gz \
    $7/flash/$pattern.extendedFrags.fastq.gz | fastq_to_fasta -o \
    $7/flash/$pattern.fa
done

# fastqc
fastqc -t $CPU -o $7/FastQC/flash $7/flash/*_flash.*.f*q.gz

# OTU picking (16S)
for f in $(find $7/flash -type f -name "*.fa"); 
do
  pick_closed_reference_otus.py -i $f \
  -r /data/genome/greengenes/fasta/97_otus.fasta \
  -o $7/OTU \
  –t /data/genome/greengenes/taxonomy/97_otu_taxonomy.txt \
  --parallel -O $CPU
done
