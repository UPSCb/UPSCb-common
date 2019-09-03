#!/bin/bash -l
#SBATCH -p core
#SBATCH -c 1
#SBATCH --mail-type=ALL

set -ex

task=$1

shift;

case "$task" in
  aggr2rmt)
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -i $1 -g $2
  ;;
  aggregate)
    gen=$1
    shift
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -g $gen -i $@
    # rename $@ to create the bin file
  ;;
  anova2el)
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -i $1 -g $2
  ;;
  ccm)
    nrow=`wc -l $1`
    ncol=`head -1 $1 | wc -w`
    ncol=`expr $ncol - 1`
    ~bastian/Git/RMTGeneNet/ccm $1 $nrow $ncol
  ;;
  el2bin)
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -i $1 -g $2 $3
  ;;
  lm2el)
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -i $1 -g $2
  ;;
  rmm)
    ~bastian/Git/RMTGeneNet/rmm -b 1 -i $1
  ;;
  threshold)
    ~bastian/Git/geneNetworkR/src/util/bin/threshold -i $1 -s 0.01 -H 0.99 -L 0.65 --trace
  ;;
  view)
    ~bastian/Git/geneNetworkR/src/util/bin/nutil --task $task -i $1 -g $2
  ;;
  ?)
    echo "No such task $task. Exiting."
    exit 1;
  ;;
esac

