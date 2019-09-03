#!/bin/bash -l

#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -p core

usage(){
    echo "runMetaxa2.sh <forward_reads> <reverse_reads> <output_directory>"
    exit l
}

frw_reads=$1
rev_reads=$2
out_dir=$3

if [ ! -e $frw_reads ]
    then
    usage
fi

if [ ! -e $rev_reads ]
    then
    usage
fi

if [ ! -d $out_dir ]
    then
    usage
fi

cd $out_dir

module load bioinfo-tools metaxa2

DBDIR=$(dirname $(which metaxa2))/metaxa2_db

metaxa2 -g ssu -1 $frw_reads -2 $rev_reads -o $out_dir -d $DBDIR --cpu 2
