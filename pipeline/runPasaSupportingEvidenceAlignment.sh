#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 16
#SBATCH --mem=100G
#SBATCH -t 4-00:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools pasa

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <config file> <genome fasta> <evidence fasta file> <out dir>

    Options:
             -c the number of CPU to use (default to 16)
             -e the step at which to end the pipeline (that step is NOT run) 
             -f the name of the file containing FL accessions
             -s the step at which to restart the pipeline
             -r DROP THE DATABASE, USE ONLY FOR CLEAN RESTART
             -I the maximal intron length size (default to 11000)
    
    Note: if -s is not provided the -C option is used instead, i.e. a new database is created
          if -e is set, no preprocessing step is done, i.e. no directory setup and no seqclean
" 
    exit 1
}

OPTIONS=
CPU=16
FLIDs=
RESET=
INTRONSIZE=11000

## get the options
while getopts c:e:f:s:rI: option
do
    case "$option" in
	c) CPU=$OPTARG;;
	e) OPTIONS="-e $OPTARG $OPTIONS";;
	f) FLIDs=$OPTARG;;
	r) RESET="reset";;
  s) OPTIONS="-s $OPTARG $OPTIONS";;
  I) OPTIONS="-I $OPTARG $OPTIONS";;
  \?) ## unknown flag
    usage;;
  esac
done
shift `expr $OPTIND - 1`

## Check if any option has been provided or if a reset has been triggered
if [[ ! -z $RESET || -z $OPTIONS ]]; then
    OPTIONS="-C"
fi 

## we get five files and a dir as argument
if [ $# != 4 ]; then
    echo "This function takes the config file, the genome fasta file, the evidence fasta file and an output dir as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing config file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing fasta file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing fasta file"
    usage
fi

if [ ! -d $4 ]; then
    echo "The fourth argument needs to be an existing directory"
    usage
fi

## cd to the output dir
cd $4

if [ "$OPTIONS" == "-C" ]; then
    ## link the files
    ln -sf $1 config
    ln -sf $2 genome
    ln -sf $3 transcripts.fasta
    
    ## run seqclean
    seqclean transcripts.fasta -c $CPU
fi    

## if we don't have accessions - assume that we are given FL sequences only
if [ -z $FLIDs ]; then
    ## extract the accessions
    $PASAHOME/misc_utilities/accession_extractor.pl < transcripts.fasta > FL_accs.txt
else
  ln -sf $FLIDs FL_accs.txt
fi

if [ ! -z $RESET ]; then
    OPTIONS="$OPTIONS -r"
fi

## run PASA
$PASAHOME/scripts/Launch_PASA_pipeline.pl -c config $OPTIONS -R -g genome -t transcripts.fasta.clean --ALIGNERS blat,gmap -u transcripts.fasta -T -f FL_accs.txt --TRANSDECODER --CPU $CPU
