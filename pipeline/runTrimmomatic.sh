#!/bin/bash
#SBATCH -p node 
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

## abort on error
set -ex

## TODO add an option for no clipping!

## usage
usage(){
echo >&2 \
"Usage: 
    Paired end: $0 <fwd fastq file> <rev fastq file> <output dir> [trimming options]
    Single end: $0 -s <fastq file> <output dir> [trimming options]

Options:
    -c      clipping file and settings
    -p      number of threads to use
    -q      use illumina quality (+64 offset), default to sanger now (+33 offset)!
    -s      single end reads
    -t      add a trim log (defaut no trimlog anymore)
    -v      verbose output


Trimming options:
    Trimming defaults to 'SLIDINGWINDOW:5:20 MINLEN:50'
    If you change the default, you need to provide the COMPLETE trimming option again!!!
    e.g. to use a 30 quality threshold for the sliding window, provide: SLIDINGWINDOW:5:30 MINLEN:50.
    Clipping defaults to 'ILLUMINACLIP:\"$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa\":2:30:10:2:keepBothReads'

"
    exit 1
}

## check env var
if [ -z $TRIMMOMATIC_HOME ]; then
    echo "The TRIMMOMATIC_HOME environment variable needs to be set. Load the T(t)rimmomatic module"
    usage
fi

## options
clip=
single_end=0
phred="-phred33"
thread=20
trimlog=0
verbose=0
while getopts "c:p:qstv" opt; do
    case $opt in
	c) clip=$OPTARG;;
	p) thread=$OPTARG;;
	q) phred="-phred64";;
        s) single_end=1;;
	t) trimlog=1;;
	v) verbose=1;;
        \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

if [ $verbose -eq 1 ]; then
    echo "Options are to use $thread CPUs"
    if [ $single_end -eq 0 ]; then
	echo "for paired-end trimming."
    else
	echo "for single-end trimming."
    fi
fi

## check the arguments
if [ $single_end -eq 0 -a $# -lt 3 ] || [ $single_end -eq 1 -a $# -lt 2 ]; then
    echo "Given the provided option, the number of argument is incorrect."
    usage
fi

## the clip default
if [ -z $clip ]; then
    if [ $single_end -eq 0 ]; then
	clip=ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads
    else
	clip=ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10:2:keepBothReads
    fi
fi

## the trim default
trim="SLIDINGWINDOW:5:20 MINLEN:50"

## check file 1
if [ ! -f "$1" ]; then
    echo "The first argument must be the valid file name of the forward fastq file."
    usage
fi
fwd=$1
shift

if [ $verbose -eq 1 ]; then
    echo "Forward file is $fwd"
fi

## create the pattern
if [ $single_end -eq 0 ]; then
    pattern=`basename ${fwd//_1.f*q.gz//}`
else
    pattern=`basename ${fwd//.f*q.gz//}`
fi

## Paired end
if [ $single_end -eq 0 ]; then

    ## check file 2
    if [ ! -f "$1" ]; then
        echo "The second argument must be the valid file name of the reverse fastq file."
        usage
    fi
    rev=$1
    shift

    if [ $verbose -eq 1 ]; then
	echo "Reverse file is $rev"
    fi
fi

## check dir
if [ ! -d "$1" ]; then
    echo $1
    echo "The third argument must be an existing directory."
    usage
fi
out=$1
shift

if [ $verbose -eq 1 ]; then
    echo "Output dir is $out"
fi

if [ $verbose -eq 1 ]; then
    echo "Phred scale is $phred"
fi

log=
if [ $trimlog -eq 1 ]; then
    log="-trimlog $out/$pattern.log"
    echo "Trim log is $log"
fi

if [ $# -gt 0 ]; then
    trim=$@
fi

if [ $verbose -eq 1 ]; then
    echo "Trimming parameters are $trim"
fi

## PE
if [ $single_end -eq 0 ]; then
    
    ## start the job
    java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE -threads $thread $phred $log $fwd $rev $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq $clip $trim

    if [ $trimlog -eq 1 ]; then
	printf "%s\0%s\0%s\0%s\0%s" $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq $out/$pattern.log | xargs -0 -I {} -P $thread gzip -f {}
    else
	printf "%s\0%s\0%s\0%s" $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq | xargs -0 -I {} -P $thread gzip -f {}
    fi

else
    java -jar $TRIMMOMATIC_HOME/trimmomatic.jar SE -threads $thread $phred $log $fwd $out/${pattern}_trimmomatic.fq $clip $trim
    if [ $trimlog -eq 1 ]; then
	printf "%s\0%s" $out/${pattern}_trimmomatic.fq $out/$pattern.log | xargs -0 -I {} -P $thread gzip -f {}
    else
	gzip -f $out/${pattern}_trimmomatic.fq	
    fi
fi

