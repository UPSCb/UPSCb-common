#!/bin/bash -l

#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -p core

usage(){
    echo "runSortmerna.sh <forward_reads> <reverse_reads> <output_directory>"
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

module load bioinfo-tools sortmerna

merged_reads=${out_dir}/merged.fastq

merge-paired-reads.sh $frw_reads $rev_reads $merged_reads

sortmerna --ref $SORTMERNADB -a 4 --log TRUE -paired_in TRUE --reads $merged_reads > ${out_dir}/summarysortmerna.log

unmerge-paired-reads.sh $merged_reads ${out_dir}/$(echo $frw_reads | cut -d '/' -f 10) ${out_dir}/$(echo $rev_reads | cut -d '/' -f 10)
