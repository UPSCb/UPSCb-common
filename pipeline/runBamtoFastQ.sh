#!/bin/bash -l

#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-05:00:00
#SBATCH --mail-type=ALL

set -e -x

##check parameters
##argument number
if [ $# != 2 ]; then
        echo -e "\e[1;31mPlease supply three arguments for this script!\e[0m"
        echo -e "\e[1;31mArgument 1 should be the input SAM or BAM file.\e[0m"
        echo -e "\e[1;31mArguments should be the output directory.\e[0m"
        exit 1
fi

## is the UPSCb env var set
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
fi

## load picard tools
module load Picard-tools

##are the files correct?
if [ ! -f $1 ]; then
        echo -e "\e[1;31mThe input file does not exist!\e[0m"
        exit 1
fi
if [ ${1: -3} != "bam" ] ; then
        echo -e "\e[1;31mThe filetype of your input file is not supported.\e[0m"
        echo -e "\e[1;31mOnly .bam is supported (case sensitive filename extension). Your file is $1\e[0m"
        exit 1
fi
if [ ! -d $2 ]; then
        echo -e "\e[1;31mThe output directory does not exist!\e[0m"
	exit 1
fi

## clean files
## run executable with suggested memory of 2GB
## the field is 1 (to report all paired reads) 
fnam=`basename ${1//.bam/}`
java -Xmx2g -jar $PICARD_HOME/SamToFastq.jar INPUT=$1 FASTQ=$2/${fnam}_1.fq SECOND_END_FASTQ=$2/${fnam}_2.fq UNPAIRED_FASTQ=$2/${fnam}_singleton.fq

