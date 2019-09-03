#!/bin/bash -l
#SBATCH -p core
#SBATCH -c 8
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools 
module load ClustalO

# usage function
usage(){
echo >&2 \
"
	Usage: $0 <fasta file to align>

"
	exit 1
}

# check values

if [ $# != 2 ]; then
    echo "This function requires only one argument."
    usage
fi

# run the command
cd $2
clustalw $1

