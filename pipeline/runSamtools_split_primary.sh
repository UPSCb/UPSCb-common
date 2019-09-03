#!/bin/bash
#SBATCH --mail-type=all
#SBATCH -p core -n 4
#SBATCH -t 1-00:00:00
#SBATCH --mem=32GB

usage(){
  echo >&2 \
  "
  Usage $0 <bam file> <out dir>
  Note: samtools >= v1.3 is expected
  "
  exit 1
}

## check if a tool is  present and is executable
toolCheck() {
    tool=`which $1 2>/dev/null`
    if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
	echo 0
    else
	echo 1
    fi
}

## version
#MAIN=1
#MAJOR=3
#MINOR=1
#versionCheck(){
#	echo $1
#  if [ $(echo $1 | awk -F. '{$1}') -ge $MAIN ] && [ $(echo $1 | awk -F. '{$2}') -ge $MAJOR ] && [ $(echo $1 | awk -F. '{$3}') -ge $MINOR ]; then
#  	echo 0
#  else
#  	echo 1
#  fi
#}

if [ $(toolCheck samtools) -eq 1 ]; then
  echo "samtools is not available"
  usage
fi

if [ ! -f $1 ]; then
  echo "The first argument needs to be a file"
  usage
fi

if [ ! -d $2 ]; then
  echo "The second argument needs to be a dir"
  usage
fi

#if [ $(versionCheck $(samtools 2>&1 > /dev/null | grep Version | awk '{print $2}')) -eq 1 ]; then
#  echo "The expected version of samtools is at least: $MAIN.$MAJOR.$MINOR"
#  usage
#fi

# extract only alignments with flag 0 or 16 (primary alignments on both strands, because the data is not strand specific)
# split alignments in the separate files by the sample names in RG tags
cd $2
samtools view -h -b -F 1797 $1 | samtools split -f "%!.%." -
