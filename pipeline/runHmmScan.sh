#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
## -A and --mail-user set in the submit job

## stop on error
set -ex

## usage
usage(){
echo >&2 \
"
	Usage: runHmmScan.sh [options] <in.hmm> <in.fa> <out dir>

	Options:
                -e report models <= this E-value threshold in output  [1e-5]

        Note:
              The script generates 4 output files and has the --acc and --noali
              parameters set as default. Here is the detail of the 6 parameters:
                -o <f>           : direct output to file <f>, not stdout
                --tblout <f>     : save parseable table of per-sequence hits to file <s>
                --domtblout <f>  : save parseable table of per-domain hits to file <s>
                --pfamtblout <f> : save table of hits and domains to file, in Pfam format <s>
                --acc            : prefer accessions over names in output
                --noali          : do not output alignments, so output is smaller
"
	exit 1
}

## load the modules
module load bioinfo-tools
module load hmmer

## options
EVALUE="1e-5"
CPU=1

## get the options
while getopts e: option
do
    case "$option" in
	e) EVALUE=$OPTARG;;
	\?) ## unknown flag
	    usage;;
    esac
done
shift `expr $OPTIND - 1`

## check the command line
if [ $# != 3 ]; then
    echo "This function takes one hmm file, one sequence file and an output directory as arguments"
    usage
fi

## check the args
if [ ! -f $1 ]; then
    echo "The first argument needs to be an existing hmm file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be an existing fasta file"
    usage
fi
nam=`basename ${2//.f*a/}`

if [ ! -d $3 ]; then
    echo "The third argument needs to be an existing directory"
    usage
fi

## run
hmmscan -E $EVALUE -o $3/$nam.txt --tblout $3/$nam.tsv --domtblout $3/${nam}-dom.tsv --pfamtblout $3/${nam}-pfam.tsv --cpu $CPU --acc --noali $1 $2

