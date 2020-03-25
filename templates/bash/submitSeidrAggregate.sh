#!/bin/bash -l

account=CHANGEME
mail=CHANGEME

# input
base=$(realpath ../data/seidr)
out=$(realpath ../data/seidr/aggregate)

# helpers
source ../UPSCb-common/src/bash/functions.sh

# modules
module load bioinfo-tools seidr-devel

# directories
if [ ! -d $base/results ]; then
  abort "Your directory structure is unexpected"
fi

if [ ! -d $base/sf ]; then
  mkdir -p $base/sf
fi

# create links
find $base/results -name "*.sf" -exec ln -sf -t $base/sf "{}" \;

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
sbatch -A $account --mem=128GB --mail-user=$mail \
  -e $out/aggregate.err -o $out/aggregate.err \
  -J aggregate ../UPSCb-common/pipeline/runSeidrAggregate.sh $out $base/sf/*.sf 
