#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools bwa

usage () {
    echo "runBwaIndex.sh -i <Genome Fasta> -o <Output dir>"
    echo
    echo "Note: If <Output dir> does not exist, it will be created"
    exit 1
}


while getopts i:o: opt
do
    case "$opt" in
	i) INDEX=$OPTARG;;
	o) OUT=$OPTARG;;
	\?)# unknown flag
	        usage;;
    esac
done

if [ -z $INDEX ] || [ -z $OUT ]; then
    echo -e "\e[91m[ERR] One or more mandatory arguments are empty\e[39m"
    usage
fi

[[ ! -f $INDEX ]] && echo -e "\e[91m[ERR] Index FASTA not found\e[39m" && usage;
[[ ! -d $OUT ]] && mkdir -p $OUT && echo -e "\e[33m[WARN] Created $OUT\e[39m";


ln -s $INDEX $OUT
BNAM=$(basename $INDEX)
echo "[INFO] Command line: runBwaIndex.sh $@"
IVDIR=$(pwd)
echo "[INFO] Invoked from: $IVDIR"
bwa index ${OUT}/$BNAM
