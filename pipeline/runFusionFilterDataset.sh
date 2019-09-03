#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools star-fusion RepeatMasker blast bowtie

# OPTIONS
CPU=32

# usage function
usage(){
echo >&2 \
"
	Usage: $0 <genome fasta> <star index> <genome gtf> <rm engine> <rm lib> <out dir>
"
	exit 1
}

# check the arguments
if [ ! -f $1 ]; then
	echo "The genome fasta file: $1 does not exist"
	usage
fi

if [ ! -d $2 ]; then
	echo "The genome directory: $2 does not exist"
	usage
fi

if [ ! -f $2/Genome ]; then
	echo "The genome directory: $2 does not seem to be a valid STAR index directory"
	usage
fi

if [ ! -f $3 ]; then
	echo "The genome gtf file: $3 does not exist"
	usage
fi

case "$4" in
  hmmer);;
  rmblast);;
  wublast);;
  *)usage;;
esac

if [ ! -f $5 ]; then
	echo "The repeat masker library file: $5 does not exist"
	usage
fi


if [ ! -d $6 ]; then
  echo "The output directory: $6 does not exist"
  usage
fi

# run the commands
cd $6
if [ ! -f cDNA_seqs.fa ]; then
  echo "Extracting the cDNA"
	$STAR_FUSION_HOME/FusionFilter/util/gtf_file_to_cDNA_seqs.pl $3 $1 > cDNA_seqs.fa
else
  echo "Skipping cDNA extraction; file exists"
fi

if [ ! -f cDNA_seqs.fa.masked ]; then
  echo "Mask the cDNA"
  RepeatMasker -pa $CPU -e $4 -s -lib $5 -xsmall cDNA_seqs.fa
else
  echo "Skipping cDNA masking; file exists"
fi

if [ ! -f cDNA_seqs.fa.masked.nsq ]; then
  echo "Construct the blast db"
  makeblastdb -in cDNA_seqs.fa.masked -dbtype nucl
else
  echo "Skipping blast db creation; files exist"
fi

if [ ! -f blast_pairs.outfmt6 ]; then
  echo "Self-self blast"
  blastn -query cDNA_seqs.fa -db cDNA_seqs.fa.masked \
      -max_target_seqs 10000 -outfmt 6 \
      -evalue 1e-3 -lcase_masking \
      -num_threads $CPU \
      -word_size 11  >  blast_pairs.outfmt6
else
  echo "Skipping blasting ; file exists"
fi

if [ ! -f blast_pairs.gene_syms.outfmt6.gz ]; then
  echo "Renaming"
  $STAR_FUSION_HOME/FusionFilter/util/blast_outfmt6_replace_trans_id_w_gene_symbol.pl \
  cDNA_seqs.fa blast_pairs.outfmt6  | gzip > blast_pairs.gene_syms.outfmt6.gz
else
  echo "Skipping renaming ; file exists"
fi

echo "Indexing"
$STAR_FUSION_HOME/FusionFilter/prep_genome_lib.pl \
           --genome_fa $1 \
           --link_star_idx $2 \
           --gtf $3 \
           --blast_pairs blast_pairs.gene_syms.outfmt6.gz \
           --cdna_fa cDNA_seqs.fa

echo "done"
