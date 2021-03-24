#!/bin/bash -l
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
## time too for large files
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

## stop on error
set -e

## be verbose and extend the commands
set -x

## load module
module load bioinfo-tools
module load trinity

## check the options if any
MEM=40G
KMER=200
PROC=20
SINGLE=0

## usage
usage(){
echo >&2 \
"
	Usage: runDigitalNormalization.sh [options] <left fq> <right fq> <out dir> [trinity options]
	runDigitalNormalization.sh [options] -s <single fq> <out dir> [trinity options]

	Options:
                -k max kmer cov (default 200)
                -m mem requirement (default 40G) (rule of thumb 1.5G/1M read pairs)
		            -p number of threads to use (default 16)
                -s for single (SE) files

        Note:
             <left fq> and <right fq> could also be files containing a list of filenames (one per line)
"
	exit 1
}

## get the options
while getopts k:m:p:s option
do
        case "$option" in
	    k) KMER=$OPTARG;;
	    m) MEM=$OPTARG;;
	    p) PROC=$OPTARG;;
	    s) SINGLE=1;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $SINGLE -eq 0 ]; then
    if [ $# -lt 3 ]; then
	echo "This function takes two files and one directory as arguments for PE data"
	usage
    fi
else 
    if [ $# -lt 2 ]; then
	echo "This function takes one files and one directory as arguments for SE data"
	usage
    fi
fi

fwd=$1
shift
if [ ! -f $fwd ]; then
    echo "The first argument needs to be the forward (left) fastq file"
    usage
fi

if [ $SINGLE -eq 0 ]; then
    rev=$1
    shift
    if [ ! -f $rev ]; then
	echo "The second argument needs to be the reverse (right) fastq file"
	usage
    fi
else
    dir=$1
    shift
    if [ ! -d $dir ]; then
	echo "The second argument (output dir) needs to be an existing directory"
	usage
    fi
fi

if [ $SINGLE -eq 0 ]; then
    dir=$1
    shift
    if [ ! -d $dir ]; then
	echo "The third argument (output dir) needs to be an existing directory"
	usage
    fi
fi

## do we have file lists
#firstLine=`head -1 $fwd`
#LIST=""
#if [ -f $firstLine ]; then
#	LIST="_list"
#	if [ $SINGLE -ne 0 ]; then
#		echo "When using SE reads, there is no support for file lists."
#		exit 1;
#	fi
#fi

## run trinity
if [ $SINGLE -eq 0 ]; then
    $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM $MEM --left$LIST $fwd --right$LIST $rev --output $dir --CPU $PROC --max_cov $KMER --pairs_together $@
else
    $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM $MEM --single $fwd --output $dir --CPU $PROC --max_cov $KMER $@
fi

