#! /bin/bash -l
#SBATCH -p node -n 8
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL

## 3 hours would be enough for most, but some nodes are slow ?!? => 6 is safe.
## even put one day!

## stop on error
set -e

## usage
usage () {
    echo "This function takes an index directory, one index name, two directories and two files as arguments"
    echo "Usage: sbatch runOSA.sh <index dir> <index name> <out dir> <tmp dir> <forward fastq> <reverse fastq>"
}

## we get on index, one dir and two files as input
if [ $# != 6 ]; then
    usage
    exit 1
fi

## check
echo Checking

## check the index
if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
    exit 1
fi

if [ ! -f $1/ReferenceLibrary/$2.dreflib1 ]; then
    echo "The second argument needs to be an existing index"
    usage
    exit 1
fi

if [ ! -d $3 ]; then
    echo "The third argument needs to be an existing directory"
    usage
    exit 1
fi

if [ ! -d $4 ]; then
    echo "The forth argument needs to be an existing directory"
    usage
    exit 1
fi

## unzip the files
if [ ! -f $5 ]; then
    echo "The fifth argument needs to be an existing fastq.gz file"
    usage
    exit 1
fi
f1=`basename ${5//_1.fastq.gz/_OSA_1.fastq}`

if [ ! -f $6 ]; then
    echo "The sixth argument needs to be an existing fastq.gz file"
    usage
    exit 1
fi
f2=`basename ${6//_2.fastq.gz/_OSA_2_fastq}`

## decompress them
echo Gunzipping
gunzip -c $5 > $4/$f1
gunzip -c $6 > $4/$f2

## chdir to out dir
echo Setting up
cd $3

## create some more var
dir=${f1//_trimmomatic_sortmerna_OSA_1.fastq/}
nam=${f1//_1.fastq/}

## if it exists...
if [ -f $dir/$nam.bam ]; then
    rm $dir/$nam.bam
fi

## align
echo Aligning
osa --alignrna -f $nam.ini -i 220 -s 40 -n $nam -d $dir -t 8 $4/$f1 $4/$f2 $1 $2 none > $nam.txt 2> $nam.err

## index
echo Indexing
module load bioinfo-tools
module load samtools/0.1.19
## for some reason, osa does not use the file name to name the bam, but the dir name if the names contains a dot...
if [ ! -f $dir/$nam.bam ]; then
    mv $dir/${dir//.*/}.bam $dir/$nam.bam
fi

## now index
samtools index $dir/$nam.bam

## clean up
echo Cleaning up
printf "%s\0%s" $4/$f1 $4/$f2 | xargs -0 -I {} -P 2 rm {}

##
echo Done
