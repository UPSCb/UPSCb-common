#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -A facility
#SBATCH -J plaac
#SBATCH --mail-usage=END,FAIL

# sanity
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Synopsis [options] $0 <plaac singularity image> <input protein fasta file> <output directory>
  
  Options:
    -p plot (returns per residue values)

  Note: the default returns per sequence value. If you use -p, provide a small fasta file, NOT a while proteome!
"

# vars
option=""

# options
while getopts p option
do
    case "$option" in
    p) option="-p all";;
    \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# checks
[[ $# -ne 3 ]] && abort "This script expects three arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file path to a singularity container"
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"
[[ ! -f $2 ]] && abort "The second argument needs to be an existing file path to the input protein fasta file"
[[ ! -d $3 ]] && abort "The third argument needs to be an existing directory"

# run
singularity exec $1 \
plaac.sh -i $2 $option > $3
