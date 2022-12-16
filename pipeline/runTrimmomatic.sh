#!/bin/bash
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=END,FAIL

## abort on error
set -ex

## source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
export USAGETXT="
Usage: 
    Paired end: $0 [options] <singularity file> <adapter file> <fwd fastq file> <rev fastq file> <output dir> [trimmomatic trimming options]
    Single end: $0 -s <singularity file> <adapter file> <fastq file> <output dir> [trimmomatic trimming options]

Options:
    -c      clipping file and settings
    -C      delete the input after successful completion
    -p      number of threads to use
    -q      use illumina quality (+64 offset), default to sanger now (+33 offset)!
    -s      single end reads
    -t      add a trim log (defaut no trimlog anymore)

Trimming options:
    Trimming defaults to 'SLIDINGWINDOW:5:20 MINLEN:50'
    If you change the default, you need to provide the COMPLETE trimming option again!!!
    e.g. to use a 30 quality threshold for the sliding window, provide: SLIDINGWINDOW:5:30 MINLEN:50.
    Clipping defaults to 'ILLUMINACLIP:<adapter file>:2:30:10:2:TRUE' for PE and 'ILLUMINACLIP:<adapter file>:2:30:10:2' for SE

"

# defaults
clip=
log=
phred="-phred33"
single_end=0
thread=20
trim="SLIDINGWINDOW:5:20 MINLEN:50"
trimlog=0
clean=0

# options
while getopts "c:Cp:qstv" opt; do
    case $opt in
	  c) clip=$OPTARG;;
	  C) clean=1;;
	  p) thread=$OPTARG;;
	  q) phred="-phred64";;
    s) single_end=1;;
	  t) trimlog=1;;
    \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## check the arguments
if [ $single_end -eq 0 -a $# -lt 5 ] || [ $single_end -eq 1 -a $# -lt 4 ]; then
    abort "Given the provided options, the number of argument is incorrect."
fi

[[ ! -f $1 ]] && abort "The first argument must be the trimmomatic singularity container file"
img=$1
shift

## enforce singularity
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -f $1 ]] && abort "The second argument must be the trimmomatic adapter fasta file"
adpt=$1
shift

## the clip default
if [ -z $clip ]; then
    if [ $single_end -eq 0 ]; then
	    clip=ILLUMINACLIP:$adpt:2:30:10:2:TRUE
    else
	    clip=ILLUMINACLIP:$adpt:2:30:10:2
    fi
fi

[[ ! -f $1 ]] && abort "The third argument must be the forward fastq file"
fwd=$1
shift

## create the pattern
pattern=$(basename ${fwd//.f*q.gz//})
if [ $single_end -eq 0 ]; then
  [[ ! -f $1 ]] && abort "The fourth argument must be the reverse fastq file"
  rev=$1
  shift
fi

[[ ! -d $1 ]] && abort "The last argument must be the output directory"
out=$1
shift

if [ $trimlog -eq 1 ]; then
    log="-trimlog $out/${pattern}.log"
fi

if [ $# -gt 0 ]; then
    trim=$@
fi

## PE
if [ $single_end -eq 0 ]; then
  singularity exec $img trimmomatic PE -threads $thread $phred $log $fwd $rev $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq $clip $trim
  if [ $trimlog -eq 1 ]; then
	  printf "%s\0%s\0%s\0%s\0%s" $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq $out/$pattern.log | xargs -0 -I {} -P $thread gzip -f {}
  else
	  printf "%s\0%s\0%s\0%s" $out/${pattern}_trimmomatic_1.fq $out/${pattern}_unpaired_1.fq $out/${pattern}_trimmomatic_2.fq $out/${pattern}_unpaired_2.fq | xargs -0 -I {} -P $thread gzip -f {}
  fi
  [[ $clean -eq 1 ]] && [[ -f $out/${pattern}_trimmomatic_1.fq.gz ]] && [[ -f $out/${pattern}_trimmomatic_2.fq.gz ]] && rm $fwd $rev
else
  singularity exec $img trimmomatic SE -threads $thread $phred $log $fwd $out/${pattern}_trimmomatic.fq $clip $trim
  if [ $trimlog -eq 1 ]; then
	  printf "%s\0%s" $out/${pattern}_trimmomatic.fq $out/$pattern.log | xargs -0 -I {} -P $thread gzip -f {}
  else
	  gzip -f $out/${pattern}_trimmomatic.fq
  fi
  [[ $clean -eq 1 ]] && [[ -f $out/${pattern}_trimmomatic.fq.gz ]] && rm $fwd
fi

# clean exit
exit 0
