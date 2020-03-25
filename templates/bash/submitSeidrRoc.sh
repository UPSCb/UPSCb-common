#!/bin/bash -l

account=CHANGEME
email=CHANGEME

# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
pgs=CHANGEME
ngs=CHANGEME
# example: 
# pgs=$(realpath ../goldStandard/Picea-abies_KEGG-based-positive-gold-standard.tsv)
# ngs=$(realpath ../goldStandard/Picea-abies_KEGG-based-negative-gold-standard.tsv)

bb=$(realpath ../data/seidr/aggregate/aggregated.sf)
indir=$(realpath ../data/seidr/backbone)
out=$(realpath ../data/seidr/roc)

if [ ! -d $out ]; then
  mkdir -p $out
fi

bbnam=$(basename $bb)
if [ ! -h $indir/$bbnam ]; then
  ln -sf -t $indir $bb 
fi

# find the network files
for f in $(find $indir -name "*.sf"); do
  fnam=$(basename ${f/.sf/})

  # run the roc on all
  sbatch -A $account --mail-user=$email \
  -o $out/${fnam}_roc.out -e $out/${fnam}_roc.err \
  ../UPSCb-common/pipeline/runSeidrRoc.sh $f $pgs $ngs \
  $out/${fnam}_roc.tsv
done
