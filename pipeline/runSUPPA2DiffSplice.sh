#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH --time=1-00:00:00
#SBATCH --job-name="suppa2"
#SBATCH --mail-type=END,FAIL

# sanity
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Synopsis [options] $0 <suppa2 singularity container> <event file> <psi comma separated files> <exp comma separated files> <out dir>
  
  Options:
    -c perform all samples combinations 
    -m use median instead of mean 
    -s export the tpm results 
    -p consider the conditions as paired data
    -t multiple testing correction
  
  Note: all options defaults to true, set the flag to deactivate them
"

# vars
paired="-pa"
combine="-c"
multiple="-gc"
export="-s"
median="-me"

# options
while getopts cmpst option
do
    case "$option" in
    c) combine="";;
    m) median="";;
    p) paired="";;
    s) export="";;
    t) multiple="";;
    \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

OPTIONS="$combine $median $paired $export $multiple"

# checks
[[ $# -ne 5 ]] && abort "This script expects five arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file path to a singularity container"
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"
[[ ! -f $2 ]] && abort "The second argument needs to be an existing file path to an event annotation file, .ioi or .ioe"
[[ ! -d $5 ]] && abort "The fifth argument needs to be an existing directory"

# separate the input files and test them
tmp=$(mktemp)
echo $3 | xargs -d, -I {} bash -c '[[ ! -f $(realpath $0) ]] && echo $0 is not a valid file > $1' {} $tmp && echo
[[ $(wc -l $tmp|cut -f1 -d " ") -ne 0 ]] && abort "Some psi files do not exist, check $3"
echo $4 | xargs -d, -I {} bash -c '[[ ! -f $(realpath $0) ]] && echo $0 is not a valid file > $1' {} $tmp && echo
[[ $(wc -l $tmp | cut -f1 -d" ") -ne 0 ]] && abort "Some tpm files do not exist, check $4"

# run
# singularity exec $1 \
suppa.py diffSplice \
-m empirical \
-p $(echo $3 | tr , " ") \
-e $(echo $4 | tr , " ") \
-i $2 -o $5 $OPTIONS
