#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 16
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools pasa/2.0.3

## TODO
## integrate the -I option for the max intron length size

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <config file> <genome fasta> <trinity fasta file> <trinity genome guided fasta file> <cufflinks gtf file> <out dir>

    Options:
             -c the number of CPU to use (default to 16)
             -e the step at which to end the pipeline (that step is NOT run)
             -s the step at which to restart the pipeline
             -n do not setup the DATABASE
             -r DROP THE DATABASE, USE ONLY FOR CLEAN RESTART
    Note: if -s is not provided the -C option is used instead, i.e. a new database is created
          if -e is set, no preprocessing step is done, i.e. no directory setup and no seqclean
"
    exit 1
}

OPTIONS=
CPU=16
RESET=
NODB=0

## get the options
while getopts c:e:s:rn option
do
  case "$option" in
	  c) CPU=$OPTARG;;
	  e) OPTIONS="-e $OPTARG $OPTIONS";;
	  n) NODB=1;;
	  r) RESET="reset";;
    s) OPTIONS="-s $OPTARG $OPTIONS";;
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
if [ $# != 6 ]; then
    echo "This function takes the config file, the genome fasta file, trinity fasta file, trinity genome guided fasta file, cufflinks gtf file and an output dir as argument"
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

if [ ! -f $4 ]; then
    echo "The fourth argument needs to be an existing fasta file"
    usage
fi

if [ ! -f $5 ]; then
    echo "The fifth argument needs to be an existing cufflinks file"
    usage
fi

if [ ! -d $6 ]; then
    echo "The sixth argument needs to be an existing directory"
    usage
fi

# ## create a temporary directory
# tmp=`mktemp`
# echo "Writing to the $tmp dir"

## cd there
# cd $tmp

## cd to the output dir
cd $6

if [ "$OPTIONS" == "-C" ]; then
    ## link the files
    ln -sf $1 config
    ln -sf $2 genome
    ln -sf $3 trinity
    ln -sf $4 trinityGG
    ln -sf $5 cufflinks

    ## combine the fasta
    zcat trinity trinityGG > transcripts.fasta
    
    ## TODO - add copying of large genome if needed
    ## mkdir genome.gmap
    ## ln -s /mnt/picea/storage/reference/Picea-abies/v1.0/indices/GMAP/Pabies1.0/Pabies1.0.* .
    ## rename "s/Pabies1.0/genome.gmap/" Pabies1.0.*
    
    ## extract the accessions
    $PASAHOME/misc_utilities/accession_extractor.pl < trinity > tdn.accs
    
    ## run seqclean
    seqclean transcripts.fasta
    
    ## clean up the flag if no db
    if [ $NODB -eq 1 ]; then
      OPTIONS=
    fi
    
fi

if [ ! -z $RESET ]; then
    OPTIONS="$OPTIONS -r"
fi

## run PASA
$PASAHOME/scripts/Launch_PASA_pipeline.pl -c config $OPTIONS -R -g genome -t transcripts.fasta.clean --ALIGNERS gmapl --TDN tdn.accs --cufflinks_gtf cufflinks -u transcripts.fasta -T --TRANSDECODER --CPU $CPU
