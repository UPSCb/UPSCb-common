#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --mail-type=ALL

usage(){
  echo >&2 \
  "
  This script expects one argument: the directory to clean up.
  
  WARNING: this script REMOVES files considered unneeded from a trinity run.
  DO NOT USE IT FOR ANYTHING ELSE!
  "
  exit 1
}

if [ $# != 1 ]; then
  echo "ERROR: This script expect one argument: the directory to clean"
  usage
fi

if [ ! -d $1 ]; then
  echo "ERROR: The provided directory does not exist"
  usage
fi

cd $1

if [ ! -f Trinity.fasta ]; then
  echo "This does not look like a completed trinity run or a trinity directory. Aborting"
  exit 1
else
  ## remove
  rm -rf \
  chrysalis \
  both.fa* \
  bowtie.* \
  FailedCommand \
  inchworm.* \
  iworm_scaffolds* \
  jellyfish.* \
  partitioned_reads* \
  read_partitions \
  recursive_trinity* \
  scaffolding_entries.sam \
  target* \
  tmp.iworm.fa* \
  norm_for_read_set* \
  
  ## compress diginorm
  find insilico_read_normalization* -type f -name "*.fq" | xargs -P 2 -I {} gzip {}
  
  ## compress output
  gzip Trinity.fasta
fi
