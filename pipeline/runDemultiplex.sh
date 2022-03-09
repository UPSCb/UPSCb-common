#!/bin/bash -l
#SBATCH -p rbx -n 1
#SBATCH --mail-type=END,FAIL
#SBATCH -t 2:00:00

set -eu

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
SYNOPSIS $0 [options] <out> <barcode> <fwd> [rev]

OPTIONS
  -r barcode is in the read
  -e end position of the barcode in the read
  
NOTE
  -r requires -e

"

isExec demultiplex

READ=
END=

while getopts e:r option
do
    case "$option" in
    e) END="-e $OPTARG";;
    r) READ="-r";;
	  \?) usage;;
  esac
done
shift `expr $OPTIND - 1`

[[ ! -z $READ ]] && [[ -z $END ]] && abort "-r requires -e"

[[ -z $READ ]] && [[ ! -z $END ]] && abort "-e requires -r"

[[ ! -d $1 ]] && abort "The output directory does not exist"

[[ ! -f $2 ]] && abort "The barcode file does not exist"

[[ ! -f $3 ]] && abort "The forward file does not exist"

[[ $# -lt 3 ]] || [[ $# -gt 4 ]] && abort "The script expects 3 or 4 arguments"

[[ $# -eq 4 ]] && [[ ! -f $4 ]] && abort "The reverse file does not exist"

[[ $# -eq 3 ]] && demultiplex demux $READ $END -p $1 $2 $3

[[ $# -eq 4 ]] && demultiplex demux $READ $END -p $1 $2 $3 $4
