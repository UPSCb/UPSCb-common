#!/bin/bash
#SBATCH --mail-type=all

## stop on error
set -e

## be verbose and extend the commands
set -x

## Genome
GENOME=12000000
NAME=
CONTROL=1
PEMODE=
## usage
usage(){
echo >&2 \
"
	Usage: runMacs2.sh [options] <treatment file> <control file> <out dir>

  Options:
    -c: no control file
    -g: set the genome size; default to 12,000,000
    -n: set the output name
    -p: paired-end mode
"
	exit 1
}

## get the options
while getopts cg:n:p option
do
        case "$option" in
        c)CONTROL=0;;
      g) GENOME=$OPTARG;;
      n) NAME="$OPTARG";;
	p)PEMODE="-f BAMPE";;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

#if [ $# != 3 ]; then
#    echo "This function takes two files and a dir as arguments"
#    usage
#fi

treatment=$1
shift
if [ ! -f $(readlink -f $treatment) ]; then
    echo "The first argument needs to be the treatment BAM file"
    usage
fi

if [ $CONTROL -eq 1 ]; then
  control=$1
  shift
  if [ ! -f $control ]; then
    echo "The second argument needs to be the control bam file"
    usage
  fi
fi

out=$1
shift
if [ ! -d $out ]; then
  echo "The third argument needs to be the output directory"
  usage
fi

## the name
if [ -z $NAME ]; then
	if [ $CONTROL -eq 1 ]; then
	  NAME="-n $(basename ${treatment//.bam/})-$(basename ${control//.bam/})"
	else
	  NAME="-n $(basename ${treatment//.bam/})"
	fi
else
  NAME="-n $NAME"
fi

## the control flag
if [ $CONTROL -eq 1 ]; then
  control="-c $control"
else
  control=
fi

macs2 callpeak -t $treatment $control -f BAM -g $GENOME --keep-dup auto --outdir $out $NAME --nomodel --extsize 150 --broad --broad-cutoff 0.05   



#call summits
#deconvolve subpeaks within each peak called from general procedure. 
#It's highly recommended to detect adjacent binding events.
# --call-summits

#-B --SPMR

##Running MACS2 bdgcmp to generate fold-enrichment and to remove background noise from BedGraph signal files for reported peaks

#options
# -m FE: calculate fold enrichment. Other options can be <logLR> log likelihood; <subtract> subtracting noise from treatment sample.
# -p: sets pseudocount. This number will be added to 'pileup per million reads' value. \
# Not needed while generating fold enrichment track because control lambda will always >0. \
# But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.


#macs2 bdgcmp -t $treatment $control -o $out -m FE

