#!/bin/bash -l
#SBATCH -p core -n 20
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -w picea

# stop on error, be verbose and expand the commands
set -e -x

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## load modules
module load bioinfo-tools diamond

## usage
USAGETXT=\
"
	Usage: runPLncPRO.sh <fasta file> <model file> <diamonddb> <out dir>
	
	Options:
	            -t              number of threads  
"

CPU=20

## get the options
while getopts t: option
do
    case "$option" in
	t) CPU=$OPTARG;;
        \?) ## unknown flag
            usage;;
  esac
done
shift `expr $OPTIND - 1`

# Check
if [ $# -ne 4 ]; then
    echo "This function needs 4 arguments"
    usage
fi

if [ ! -f $1 ]; then
  abort "The first argument needs to be the fasta filepath"
fi

if [ ! -f $2 ]; then
  abort "The second argument needs to be the model filepath"
fi

if [ ! -f $3 ]; then
  abort "The third argument needs to be the diamond index filepath"
fi

if [ ! -d $4 ]; then
    abort "The fourth argument (output dir) needs to be an existing directory"
fi

# run PLncPRO
#cd $2

# output filename
fnam=$(basename ${1%.*})
singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/delhomme-upscb-lncrna.simg python /opt/plncpro/prediction.py \
-t $CPU -i $1 -m $2 -d $3 -o $4 -p $4/$fnam.txt
#singularity exec delhomme-upscb-lncrna.simg python /opt/plncpro/prediction.py
#singularity exec delhomme-upscb-lncrna.simg diamond help
