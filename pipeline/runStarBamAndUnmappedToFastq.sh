#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

## stop on error, be verbose and expand the commands
set -e -x

## load necessary modules
module load bioinfo-tools samtools # Picard-tools

## usage
usage(){
echo >&2 \
"
	Usage: $0 <bam file> <unmapped fwd fq.gz> <unmapped rev fq.gz> <out dir> 
"
	exit 1
}

# no 5th arg <chimeric bam file>

## check the input
# if [ $# -lt 4 ]; then
if [ $# != 4 ]; then
  echo "This script expects 4 arguments"
  usage
fi

#  at least; a 5th (Chimeric bam file) is optional

if [ ! -f $1 ]; then
  echo "The first argument needs to be an existing BAM file"
  usage
fi

if [ ! -f $2 ]; then
  echo "The second argument needs to be an existing FASTQ file"
  usage
fi

if [ ! -f $3 ]; then
  echo "The third argument needs to be an existing FASTQ file"
  usage
fi

if [ ! -d $4 ]; then
  echo "The forth argument needs to be an existing output directory"
  usage
fi

# ## extract the files
# fnam=`basename ${1//.bam/}`
# fwd=$4/${fnam}_1.fq
# rev=$4/${fnam}_2.fq
# unp=$4/${fnam}.fq
# 
# nice idea but fails as picard does not understand that aligned reads can have 
# unmapped mate...
# java -jar /mnt/picea/Modules/apps/bioinfo/picard/1.117/SamToFastq.jar INPUT=$1 \
# FASTQ=$rev SECOND_END_FASTQ=$fwd UNPAIRED_FASTQ=$unp
# 
# ## concatenate them
# grep "\1" $unp | paste - - - - | sort -k1,1 -S 6G | tr '\t' '\n' >> $fwd
# grep "\2" $unp | paste - - - - | sort -k1,1 -S 6G | tr '\t' '\n' >> $rev

## clean
# rm $unp

# so instead the solution with Picard would be a HUGE FREAKING DETOUR
# # 1. fastq to SAM
# java -jar $PICARD_HOME/FastqToSam.jar \
# FASTQ= \
# FASTQ2= \
# OUTPUT= $out/$fnam.sam
# 
# # 2. sort
# samtools view -bS $out/$fnam.san | samtools sort -n - $out/$fnam.bam
# 
# # 3. merge
# java -jar $PICARD_HOME/MergeBamAlignment.jar \
# UNMAPPED_BAM= \
# ALIGNED_BAM= \
# REFERENCE= \
# OUTPUT= 
# 
# # 4. and finally back to fastq
# java -jar /mnt/picea/Modules/apps/bioinfo/picard/1.117/SamToFastq.jar
# INPUT= \
# FASTQ=$rev \
# SECOND_END_FASTQ=$fwd

# Since I don't care or reverting sequences on the - strand, I just
# pull the first mate with samtools, the second mate with samtools
# cat that to the unmapped and sort
fnam=`basename ${2//_Unmapped_1.fq.gz/_STAR}`
fwd=$4/${fnam}_1.fq
rev=$4/${fnam}_2.fq
## to drop the unneeded header section
zcat $2 | awk '{print $1}' > $fwd
zcat $3 | awk '{print $1}' > $rev

# first (or second) read in pair and primary
samtools view -h -f 0x0040 $1 | samtools view -F 256 - | awk '{print "@"$1"/1\n"$10"\n+\n"$11}' >> $fwd
samtools view -h -f 0x0080 $1 | samtools view -F 256 - | awk '{print "@"$1"/2\n"$10"\n+\n"$11}' >> $rev

# no chimeric as it contains split reads
# not worth the time implementing how to restore these
#if [ ! -z $5 ]; then
#  samtools view -f 0x0040 $5 | awk '{print $1"/1\n"$9"\n+\n"$11}' >> $fwd
#  samtools view -f 0x0080 $5 | awk '{print $1"/2\n"$9"\n+\n"$11}' >> $rev
#fi

cat $fwd | paste - - - - | awk '{f=$1;sub(/\/1/,"",f);print f,$1,$2,$3,$4}' | sort -k1,1 -S 6G > $4/$fnam.tmp_1
cat $rev | paste - - - - | awk '{f=$1;sub(/\/2/,"",f);print f,$1,$2,$3,$4}' | sort -k1,1 -S 6G > $4/$fnam.tmp_2
join -1 1 -2 1 $4/$fnam.tmp_1 $4/$fnam.tmp_2 > $4/$fnam

cut -f 2,3,4,5 -d" " $4/$fnam | tr ' ' '\n' | gzip > $fwd.gz
cut -f 6,7,8,9 -d" " $4/$fnam | tr ' ' '\n' | gzip > $rev.gz

rm $fwd $rev $4/$fnam $4/$fnam.tmp_1 $4/$fnam.tmp_2
