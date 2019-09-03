#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools pasa

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <config file> <gff3 file> <pasa dir> 
    Options:
             -r <database name> start the comparison (do not load nor validate the gff3) on the specified database
             -v skip the gff validation. This might be necessary for re-run (split genes keep their original IDs!) 
    Note: 
             Using the -r option IMPLIES that the gff annotation have been loaded previously!
" 
# multithreading failed (hangs)

#             -c the number of CPU to use (default to 16)    
    exit 1
}

## GLOBAL VARS
CPU="--CPU 1"
RESTART=
VALIDATE=1

## get the options
# while getopts c:r option
 while getopts r:v option
 do
     case "$option" in
# multithreading failed (hangs)
#   c) CPU="--CPU $OPTARG";;
	 r) RESTART=$OPTARG;;
	 v) VALIDATE=0;;
   \?) ## unknown flag
    usage;;
   esac
 done
 shift `expr $OPTIND - 1`

## we get two file and a dir as input
if [ $# != 3 ]; then
    echo "This function takes a config file, a gff3 file and an existing pasa dir as argument"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing config file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing gff3 file"
    usage
fi

if [ ! -d $3 ]; then
    echo "The third argument needs to be an existing directory"
    usage
    if [[ ! -f genome || ! -f transcripts.fasta.clean ]]; then
	echo "This does not look like a pasa results directory"
	usage
    fi
fi

## cd to the output dir
cd $3

if [ -z $RESTART ]; then

    if [ $VALIDATE -eq 1 ]; then
	## validate the gff
	$PASAHOME/misc_utilities/pasa_gff3_validator.pl $2
    fi

    ## run PASA
    $PASAHOME/scripts/Launch_PASA_pipeline.pl -c $1 -A -g genome -t transcripts.fasta.clean -L --annots_gff3 $2 $CPU
else
    $PASAHOME/scripts/cDNA_annotation_comparer.dbi -G genome -M $RESTART $CPU
fi
