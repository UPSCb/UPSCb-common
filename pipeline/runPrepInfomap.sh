#!/bin/bash

inf=$1
prep="/mnt/picea/home/bastian/tmp/prep"
back="/mnt/picea/home/bastian/tmp/back"
info="/mnt/picea/home/bastian/tmp/Infomap"
outdir=$2
bn=$(basename $inf)

$prep $inf map > $outdir/$bn.map
$prep $inf print > $outdir/$bn.infomap.txt
$info -z -i link-list --markov-time 0.01 $outdir/$bn.infomap.txt $outdir
$back $outdir/$bn.map $outdir/$bn.infomap.tree > $outdir/$bn.final.txt

nc=$(head -n 1 $outdir/$bn.final.txt | wc -w)
nc=$(expr $nc - 3)

>$outdir/$bn.final.h.txt
for f in $(seq 1 $nc); do echo -ne "C$f\t" >>  $outdir/$bn.final.h.txt; done
echo -e "Score\tIID\tID" >>  $outdir/$bn.final.h.txt
cat  $outdir/$bn.final.txt >>  $outdir/$bn.final.h.txt
