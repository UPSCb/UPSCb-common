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
  Synopsis [options] $0 <suppa2 singularity container> <gtf file> <out dir>
  
  Options:
    -f the format to export (ioe, default, or ioi)
    -e disable the generation of SE
    -i disable the generation of RI
    -l disable the generation of FL
    -s disable the generation of SS
    -x disable the generation of MX
"

# arguments
SE=SE
SS=SS
MX=MX
RI=RI
FL=FL
format=ioe

# options
while getopts ef:ilsx option
do
    case "$option" in
    e) SE="";;
    f) format=$OPTARG;;
    i) RI="";;
    l) FL="";;
    s) SS="";;
    x) MX="";;
    \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# checks
[[ $# -ne 3 ]] && abort "This script expects three arguments"
[[ ! -f $1 ]] && abort "The first argument needs to be an existing file path to a singularity container"
[[ -z ${SINGULARITY_BINDPATH:-} ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"
[[ ! -f $2 ]] && abort "The second argument needs to be an existing file path to a gtf annotation file"
[[ ! -d $3 ]] && abort "The third argument needs to be an existing directory"

# run
singularity exec $1 \
suppa.py generateEvents \
-i $2 \
-o $3/$(basename $3) \
-e $SE $SS $MX $RI $FL \
-f $format
