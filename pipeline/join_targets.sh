#!/bin/bash

dir=/mnt/picea/home/katja/aspen_sRNA/Potra/ShortStack/targets

cat $dir/psRNATarget_*_Potra_4.txt >> $dir/targets_all_Potra_4.txt
grep "miRNA_[0-9]" $dir/targets_all_Potra_4.txt > $dir/targets_all_Potra.tmp
grep "miRNA_Acc." $dir/psRNATarget_21nt_Potra_4.txt > $dir/targets_all_Potra_4.txt
cat $dir/targets_all_Potra.tmp >> $dir/targets_all_Potra_4.txt
rm $dir/targets_all_Potra.tmp

cat $dir/psRNATarget_*_Potri_4.txt >> $dir/targets_all_Potri_4.txt
grep "miRNA_[0-9]" $dir/targets_all_Potri_4.txt > $dir/targets_all_Potri.tmp
grep "miRNA_Acc." $dir/psRNATarget_21nt_Potri_4.txt > $dir/targets_all_Potri_4.txt
cat $dir/targets_all_Potri.tmp >> $dir/targets_all_Potri_4.txt
rm $dir/targets_all_Potri.tmp
