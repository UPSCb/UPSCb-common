#!/bin/bash -l

#SBATCH -p node
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH -n 8

##be verbose and stop on error
set -ex

##usage
usage(){
    echo >&2 \
"Usage: $0 [options] <input sample fasta> <output prefix> <centroid database fasta>

A database of globally clustered OTUs needs to be available, create BEFORE running this script by running runVsearchDatabase.sh

Options:
      -t  type of input data (16S or ITS, default: 16S) 
      -i  sequence identity (if not specified, it will be set to the default values of 0.97 for 16S and 0.95 for ITS)
      -q  write non-matching query sequences to separate fasta file

Note:
  You need to set the UPSCb env. variable to your UPSCb git checkout directory
"
    exit 1
}

##load modules
module load bioinfo-tools vsearch Qiime/1.9.0 blast/2.2.26

##defaults
typ=16S
ident=0.97

##options
query=0
while getopts t:i:q opt;
do
    case $opt in
	t) typ=$OPTARG;;
	i) ident=$OPTARG;;
	q) query=1;;
	\?) usage;;
    esac
done

shift `expr $OPTIND - 1`

##create variables depending on input options
if [ $typ == "ITS" ]; then
    if [ $ident == "0.97" ]; then
	ident=0.95
    fi
fi

##arguments
if [ $# != 3 ]; then
    echo "This function takes 3 arguments, the input fasta file, the output directory and a centroid (nonchimeras.fa) database fasta file"
    usage
fi

#if [ ! -f $1 ]; then
#    echo "The first argument must be a valid fasta file"
#    usage
#fi

if [ ! -f $3 ]; then
    echo "The third argument must be a valid fasta file"
    usage
fi


if [ ! -d $2 ]; then
    mkdir -p $2
fi

##variable to create ouput filename containing sampleID
BN=$(basename $1)
ID=${BN/demultiplex_/dm_}
FFF=${ID/_trimmomatic_flash.fa/_}


##commands

###########OTU POPULATION######################################
#all reads in the sample from the input file are searched against the centroid file containing ALL OTUs from the whole library/libraries, determined beforehand with runVsearchDatabase.sh

if [ $query -eq 0 ]; then
   vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --usearch_global $1 --db $3 --id $ident --strand both --otutabout $2/${FFF}OTUtable.txt --dbmatched $2/${FFF}db_out.fa --matched $2/${FFF}query_out.fa --fastapairs $2/${FFF}pairs.fa --sizeout
   else
       vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --usearch_global $1 --db $3 --id $ident --strand both --otutabout $2/${FFF}OTUtable.txt --dbmatched $2/${FFF}FASTAout.fa --matched $2/${FFF}query_out.fa --fastapairs $2/${FFF}pairs.fa --notmatched $2/${FFF}notmatched.fa --sizeout
fi
   
