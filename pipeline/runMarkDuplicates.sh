#!/bin/bash -l

#SBATCH -t 1-00:00:00
#SBATCH -p core
#SBATCH -n 1

set -e

module load java
module load bioinfo-tools
module load samtools
#module load picard
module load Picard-tools

THREADS=1

if [ -z $PICARD_HOME ]; then
    echo >&2 "Could not find picard tools"
    exit 1
fi

# default
JavaMem=6G

if [ $# -ne 2 ]; then
    echo "Usage: $0 <BAM file> <output directory>" 1>&2
    exit 1
fi

if [ ! -f $1 ]; then
    echo "Could not find BAM file '$1'" 1>&2
    exit 1
fi
inbam=$1

if [ ! -d $2 ]; then
    echo "Could not find directory '$2'" 1>&2
    exit 1
fi
outdir=$2

sname=`basename "${inbam/_[st]*[st]*_STAR.bam/}"`
name_out=`basename "${inbam/.bam/}"`

# Run MarkDuplicates
java -Xmx${JavaMem} -XX:ParallelGCThreads=$THREADS -jar $PICARD_TOOLS_DIR/picard.jar MarkDuplicates \
    ASSUME_SORTED=true \
    INPUT=$inbam \
    OUTPUT=$outdir/${name_out}_mkdup.bam \
    METRICS_FILE=$outdir/${sname}_mkdup.metrics \
    CREATE_INDEX=true

