#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL

usage(){
echo >&2 \
"
	Usage: $0 <genome>
	
	Arguments:
                genome: The genome STAR index directoy
"
	exit 1
}

## load module
module load bioinfo-tools star

## check that the genome exists
if [ $# != 1 ]; then
    echo "This function takes one argument the STAR genome"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to point to a valid STAR index directory"
    usage
fi

## load the genome
STAR --genomeDir $1 --genomeLoad LoadAndExit

