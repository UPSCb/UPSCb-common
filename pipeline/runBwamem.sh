#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools bwa/0.7.17 samtools

set -ex

usage () {
    echo >&2 "
    Usage: runBwamem.sh -f <Forward read> -r <Reverse read> -i <Index Fasta> -o <Output dir>
    
    Optional arguments:
    -t <INT>     ----- Number of mapping and sorting threads
    -I           ----- Index the fasta file if no index is present
    -O <string>  ----- optional arguments to pass to bwa mem (must be quoted)
    -s           ----- single end alignment
    "                                                
    echo -e "\e[91mNOTE: You need to reserve memory for -t <INT> * 768M in SLURM\e[39m"
    exit 1
}


OPTS=""
SINGLE=0

while getopts f:i:r:t:o:O:Is opt
do
    case "$opt" in
	f) FW=$OPTARG;;
	i) INDEX=$OPTARG;;
	r) RV=$OPTARG;;
	t) THREADS=$OPTARG;;
	o) OUT=$OPTARG;;
	I) DOINDEX=1;;
	O) OPTS=$OPTARG;;
	s) SINGLE=1;;
	\?)# unknown flag
	        usage;;
    esac
done

if [ $SINGLE -eq 1 ]; then
  RV="dummy"
fi

if [ -z $FW ] || [ -z $RV ] || [ -z $INDEX ] || [ -z $OUT ]; then
    echo -e "\e[91m[ERR] One or more mandatory arguments are empty\e[39m"
    usage
fi

#echo $PPID

if [ ! -f $FW ]; then
  echo -e "\e[91m[ERR] Forward FASTQ not found\e[39m"
  usage
fi
if [ $RV != "dummy" ] && [ ! -f $RV ]; then
  echo -e "\e[91m[ERR] Reverse FASTQ not found\e[39m"
  usage
fi
if [ ! -f $INDEX ]; then
  echo -e "\e[91m[ERR] Index FASTA not found\e[39m"
  usage
fi

if [ ! -d $OUT ]; then 
  mkdir -p $OUT
  echo -e "\e[33m[WARN] Created $OUT\e[39m"
fi

if [ $DOINDEX ]; then
    NEEDS_INDEX=0
    for suffix in .ann .pac .bwt .sa .amb .fai; do
	if [ ! -f ${INDEX}$suffix ]; then
	    NEEDS_INDEX=1
	fi
    done

    if [ NEEDS_INDEX ]; then
	bwa index $INDEX
    fi
fi

[[ -z $THREADS ]] && THREADS=1;

END="_1.f*q"
if [ $RV == "dummy" ]; then
  END=".f*q"
  RV=""
fi

NAM=$(basename $FW)
if [[ $NAM =~ \.gz$ ]]; then END="${END}.gz"; fi

bwa mem $OPTS -t $THREADS $INDEX $FW $RV | samtools view -bT $INDEX - | samtools sort -@ $THREADS -o ${OUT}/${NAM/$END/.sorted.bam}

samtools index ${OUT}/${NAM/$END/.sorted.bam}
