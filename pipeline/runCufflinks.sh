#! /bin/bash -l
#SBATCH -p node -n 20
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
##
echo Loading

## stop on error
## and be verbose
set -ex

## load modules
module load bioinfo-tools cufflinks

## define a function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <out dir> <genome fasta> <alignment bam> <gtf file>

    Options:
             -c the number of CPU to use (default to 20)
             -i the max intron length (default to 11000)
             -l the library type (see 'cufflinks -h' for options, default to 'fr-unstranded')
             -n do not use --frag-bias-correct (which optimise the count estimates 
                but takes longer and requires the genome). This is good if cufflinks is run for
                assembly (e.g. for PASA)
" 
    exit 1
}

CPU="-p 20"
LIBTYPE="--library-type fr-unstranded"
INTRON="--max-intron-length 11000"
ESTIMATE=1

## get the options
while getopts c:i:l:n option
do
    case "$option" in
	c) CPU="--p $OPTARG";;
	i) INTRON="--max-intron-length $OPTARG";;
	l) LIBTYPE="--library-type $OPTARG";;
	n) ESTIMATE=0;;
        \?) ## unknown flag
            usage;;
  esac
done
shift `expr $OPTIND - 1`

##
echo Checking

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function takes a directory and three files as arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing file"
    usage
fi

if [ ! -f $4 ]; then
    echo "The forth argument needs to be an existing file"
    usage
fi

##
echo Setting up

## go to the output dir
cd $1

## get the name
nameIn=${3%.bam}
nameOut=${nameIn##*/}

OUT_PATH="$1/$nameOut"
if [[ ! -d $OUT_PATH ]]; then
	mkdir $OUT_PATH
fi

# prep options
FBC=
if [ $ESTIMATE -eq 1 ]; then
  FBC="--frag-bias-correct $2"
fi

##
echo Starting

## -v stands for verbose without the progress bar
cufflinks -v --3-overhang-tolerance 4000 --no-update-check \
--multi-read-correct $FBC $LIBTYPE $INTRON $CPU -o $OUT_PATH -g $4 $3 
#> $1/$nameOut-cufflinks.out 2> $1/$nameOut-cufflinks.err

## not needed anymore
#mv $OUT_PATH/genes.fpkm_tracking $OUT_PATH/$nameOut\_genes.fpkm_tracking
#mv $OUT_PATH/isoforms.fpkm_tracking $OUT_PATH/$nameOut\_isoforms.fpkm_tracking
#mv $OUT_PATH/skipped.gtf $OUT_PATH/$nameOut\_skipped.gtf
#mv $OUT_PATH/transcripts.gtf $OUT_PATH/$nameOut\_transcripts.gtf

##
echo Done

