#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools InfoMap

# usage function
usage(){
echo >&2 \
"
	Usage: $0 <pajek (.net) graph file> <out dir>

	Notes:
		The script only accept parjek formatted graph files
"
	exit 1
}

# check the arguments
if [ ! -f $1 ]; then
	echo "The pajek graph file: $1 does not exist"
	usage
fi


if [ "${1##*.}" != ".net"]; then
  echo "The graph file needs to be a pajek file, i.e. have a .net extension"
  usage
fi


if [ ! -d $2 ]; then
  echo "The output directory: $2 does not exist"
  usage
fi

# run the command
cd $2
Infomap -i pajek -u $1 .
