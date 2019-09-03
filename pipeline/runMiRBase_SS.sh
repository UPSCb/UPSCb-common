#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools bwa/0.7.10

dir_ref=/mnt/picea/storage/reference/miRBase/v21/indices/BWA
dir_seq=/mnt/picea/home/katja/aspen_sRNA/Potra/ShortStack/miRBase
outdir=/mnt/picea/home/katja/aspen_sRNA/Potra/ShortStack/miRBase
mkdir -p $outdir

bwa aln $dir_ref/hairpin_T.fa $dir_seq/*.precursor_T.fa > $outdir/Potra_SS_hairpin.sai
bwa samse -n 200 $dir_ref/hairpin_T.fa $outdir/Potra_SS_hairpin.sai $dir_seq/*.precursor_T.fa > $outdir/Potra_SS_hairpin.sam

bwa aln $dir_ref/hairpin_T.fa $dir_seq/*.mature_T.fa > $outdir/Potra_SS_miRNA.sai
bwa samse -n 200 $dir_ref/hairpin_T.fa $outdir/Potra_SS_miRNA.sai $dir_seq/*.mature_T.fa > $outdir/Potra_SS_miRNA.sam

grep -v "@" $outdir/Potra_SS_hairpin.sam | awk '{ if($3 != "\*") print $0 }' > $outdir/Potra_SS_hairpin.mapped
grep -v "@" $outdir/Potra_SS_miRNA.sam | awk '{ if($3 != "\*") print $0 }' > $outdir/Potra_SS_miRNA.mapped
