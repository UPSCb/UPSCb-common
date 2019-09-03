#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 0-05:00:00
#SBATCH --mail-type=AL

####
#	Post filtering for the long non coding RNA pipeline
####

usage="
	The script will filter for full length that are at least 40% of their mRNA, to be consistent with PASA. Further filterint steps are done in R
	
	runPythonLncRNAFiltering.py <frameDP folder> <pasa fasta> <output>

	1. First merge all output files from frameDP
	2. Then add all the changes to a new merged PASA_and_frameDP.gff3 from pasa dump gff and the frameDP result
	3. Check for full length status and >40% CDS
	4. Continue in R

"

### revieve arguments

## check
if [ $# -gt 3 ] || [ $# -lt 3 ]; then
	echo "Not correct number of arguments
	"
	echo $#
    usage
fi

if [ ! -f $1 ] ; then
    echo "The first argument should be a frameDP folder." $1
    usage
fi

if [ ! -d $2 ] ; then
    echo "The second argument should be the pasa_dump.gff3" $2
    usage
fi

## arguments

fdp_folder=$1
pasa_fasta=$2
output=$3

## merge all temp_gene files from the frameDP folder

cat $fdp_folder/workdir/FDP*/temp_gene*/*.gff3 > $fdp_folder/merged_frameDP.gff3

## This script will generate a new fasta file as well as the gff3 file with the new annotations, the fasta file can be used to run additional BLAST.
python $UPSCb/src/python/novel_genes/merge_fa_and_frameDP_gff.py $pasa_fasta $fdp_folder/merge_frameDP.gff3 --gff3 $output".gff3"

## check for full length status
python $UPSCb/src/python/novel_genes/check_fl_status.py -pc 0.4 -fl True -g $fdp_folder/fl_genes.txt

echo "Continue you analysis in R by using the fl_genes.txt file"