#!/bin/bash -l
#SBATCH --mail-type=all

usage () {
echo "Usage:"
echo "runRmarkdown.sh <my_R_script.R>"
echo
}

if [ ! $# == 1 -o ! -f $1 ]; then
	usage
	exit 1
fi

module load R
Rscript -e "library(rmarkdown); render(commandArgs(TRUE))" $1
