#!/bin/bash
#SBATCH --mail-type=all
#SBATCH -p node -n 16
#SBATCH --mem=128G
#SBATCH -t 4-00:00:00

# modules
#module load bioinfo-tools ShortStack

# defaults
CPU=16
MEM=120G
FORMAT="--readfile"
# usage
usage(){
  echo >&2 \
  "
    Usage: $0 [options] outdir genome file(s)
    
    Note: multiple files can be provided space separated
    
    Options:
      -b input is a bamfile
      -c the number of CPUs default to 16
      -m the memory for sorting default to 120G
  "
  exit 1;
}

# process the options
while getopts "bc:m:" opt; do
    case $opt in
      b) FORMAT="--bamfile";;
	    c) CPU=$OPTARG;;
	    m) MEM=$OPTARG;;
      \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

# check arguments
if [ $# -lt 3 ]; then
  echo "This function expects at least 3 arguments"
  usage
fi

if [ ! -d $1 ]; then
  echo "The first argument needs to be an existing directory"
  usage
fi
out=$1/`date +%Y%m%d`

if [ ! -f $2 ]; then
  echo "The second argument needs to be an existing bowtie index"
  usage
fi
bwt=$2

shift 2

for file in "$@"; do
  if [ ! -f $file ]; then
    echo "The read file $file does not exists"
    usage
  fi
done

# run
ShortStack --bowtie_cores $CPU --sort_mem $MEM --dicermin 18 \
--mismatches 0 --nostitch --outdir $out --genomefile $bwt $FORMAT $@ 
