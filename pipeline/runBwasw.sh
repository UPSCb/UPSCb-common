#!/bin/bash -l
#SBATCH --mail-type=all

module load bioinfo-tools bwa/0.7.12 samtools

usage () {
    echo "runBwasw.sh -f <Forward read> -r <Reverse read> -i <Index Fasta> -o <Output dir>"
    echo
    echo "Optional arguments:"
    echo "                       -t <INT>     ----- Number of mapping and sorting threads"
    echo "                       -I           ----- Index the fasta file if no index is present"
    echo "                       -O <string>  ----- optional arguments to pass to bwa mem
                                                    (must be quoted)"
    echo
    echo -e "\e[91mNOTE: You need to reserve memory for -t <INT> * 768M in SLURM\e[39m"
    exit 1
}


OPTS=""

while getopts f:i:r:t:o:O:I opt
do
    case "$opt" in
	f) FW=$OPTARG;;
	i) INDEX=$OPTARG;;
	r) RV=$OPTARG;;
	t) THREADS=$OPTARG;;
	o) OUT=$OPTARG;;
	I) DOINDEX=1;;
	O) OPTS=$OPTARG;;
	\?)# unknown flag
	        usage;;
    esac
done

if [ -z $FW ] || [ -z $RV ] || [ -z $INDEX ] || [ -z $OUT ]; then
    echo -e "\e[91m[ERR] One or more mandatory arguments are empty\e[39m"
    usage
fi

echo $PPID

[[ ! -f $FW ]] && echo -e "\e[91m[ERR] Forward FASTQ not found\e[39m" && usage;
[[ ! -f $RV ]] && echo -e "\e[91m[ERR] Reverse FASTQ not found\e[39m" && usage;
[[ ! -f $INDEX ]] && echo -e "\e[91m[ERR] Index FASTA not found\e[39m" && usage;
[[ ! -d $OUT ]] && mkdir -p $OUT && echo -e "\e[33m[WARN] Created $OUT\e[39m";

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
NAM=$(basename $FW)
if [[ $NAM =~ \.gz$ ]]; then END="${END}.gz"; fi


bwa bwasw $OPTS -t $THREADS $INDEX $FW $RV | samtools view -bT $INDEX - | samtools sort -@ $THREADS - -f ${OUT}/${NAM/$END/.sorted.bam}


