#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL

## stop on error
set -ex

## modules
module load bioinfo-tools pasa

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 <config file> <out dir> 

    Alt: Add cufflinks as a third argument to add a cufflinks gtf
    Note: The out dir should be an existing PASA results directory with an existing transcripts.fasta.clean file
          And if the cufflinks argument is added a link or file cufflinks must exist in the PASA directory.
" 
    exit 1
}

## we get one file and one dir as input
if [ $# != 2 | $# != 3 ]; then
    echo "This function takes one config file and one out dir"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing config file"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory file"
    usage
fi

if [ $# == 3 ] && [ $3 == "cufflinks" ]; then
    cufflinks="--cufflinks_gtf cufflinks"
elif [ $# == 2 ]; then
    echo ""
else
    echo "The third argument needs to be 'cufflinks'"
    usage
fi

## go into the directory
cd $2

## check if the file exists
if [ ! -f transcripts.fasta.clean ]; then
    echo "This does not look like a PASA results directory."
    usage
fi

## execute
$PASAHOME/scripts/build_comprehensive_transcriptome.dbi -c $1 -t transcripts.fasta.clean $cufflinks
