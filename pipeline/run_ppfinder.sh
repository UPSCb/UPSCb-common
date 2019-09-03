#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools ppfinder

# usage function
usage(){
echo >&2 \
"
	Usage: $0 
		Parameters can be changed in parameter.file
"
	exit 1
}

# run the command
cd $2
ppfinder


