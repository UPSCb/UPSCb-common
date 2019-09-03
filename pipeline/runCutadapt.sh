#!/bin/bash -l

set -e

usage(){
echo >&2 \
"
Usage: runCutadapt.sh [option} <adapter> <infile> <out dir>

Options:
        -a ADAPTER    :    This adapter file or sequence will be
                           considered ligated to the 3' end.
        -b ADAPTER    :    This adapter file or sequence will be
                           considered ligated to both ends.
        -g ADAPTER    :    This adapter file or sequence will be
                           considered ligated to the 5' end.
"
exit 1
}

if [ $# -lt 3 ] || [ $# -gt 5 ];then
    echo "Enter at least 3, but no more than 5 parameters."
    usage
    exit 1
fi

IN=0
appends=""
while getopts a:b:g: option
do
    case "$option" in
	a)
	   if [ ! -f $OPTARG ];then
	       echo "Trimming file was not found or is not a file. exiting."
	       exit 1
	   fi
	   appends="-a $OPTARG"
	   IN=$((IN + 1));;
	b)
	   if [ ! -f $OPTARG ];then
	       echo "Trimming file was not found or is not a file. exiting."
	       exit 1
	   fi
	   appends="-b $OPTARG $appends"
	   IN=$((IN + 1));;
	g)
	   if [ ! -f $OPTARG ];then
	       echo "Trimming file was not found or is not a file. exiting."
	       exit 1
	   fi
	   appends="-g $OPTARG $appends"
	   IN=$((IN + 1));;
	\?)
	    usage;;
    esac
done
IN=$((IN * 2))
shift $IN

#get other args
infile=$1
outdir=$2

#Checks for binaries
if [ -z `which cutadapt` ]; then
    echo "No local install of cutadapt detected. Trying modules"
    module load bioinfo-tools cutadapt
    if [ -z `which cutadapt` ]; then
	echo "Still no cutadapt detected. exiting"
	exit 1
    fi
    echo "Module loaded. You are using"
    cutadapt --version
fi

#Checks for files
if [ ! -f $input ];then
    echo "Input file not found or not a file. exiting."
    exit 1
fi

if [ ! -d $outdir ];then
    echo "Output directory not found or not a directory. Making it"
    mkdir -p $outdir #Nico takes the blame for all harm this causes :P
fi

#call
nam=`basename $infile`
out=$(sed -e 's/.fq/_cut.fq.gz/g' -e 's/.fq.gz/_cut.fq.gz/g' -e 's/.fastq.gz/_cut.fq.gz/g' -e 's/.fastq/_cut.fq.gz/g' <<< "$nam")
#remove double whitespace and spaces in call
call="$appends $infile -o $outdir$out"
call=${call//\ \/}
call=${call//\/\//\/}
cutadapt $call

