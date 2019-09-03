#!/bin/bash -l
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 8
## time too for large files
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

## stop on error
set -e

## be verbose and extend the commands
set -x

## check the options if any
KEEP=0
UNPAIRED=0
PROC=8

## local run
if [ -z $SLURM_SUBMIT_DIR ]; then
    SLURM_SUBMIT_DIR="."
fi

## usage
usage(){
echo >&2 \
"
	Usage: runSortmernaMtSSU.sh [option] <out dir> <tmp dir> <forward fastq.gz> <reverse fastq.gz>
	
	Options:
                -m keep the mRNA
		-p number of threads to be used (default 8)
                -u single end data (in that case only the forward fastq is needed)

        Details: $0 extracts the mitochondrial small subunit (mtSSU) reads from a fastq file whether paired
                 or not (see the -u option)
"
	exit 1
}

## get the options
while getopts mp:u option
do
        case "$option" in
	    m) KEEP=1;;
	    p) PROC=$OPTARG;;
	    u) UNPAIRED=1;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

##
echo Setting up

## set some env var
## this location is not in Git anymore!
## it has to be downloaded by the user
## check the ethylene-insensitive project submitter to see
## how to set that up
if [ -z $SORTMERNADIR ]; then
    export SORTMERNADIR=$SLURM_SUBMIT_DIR/../../../data/sortmerna
fi

## set the dbs
mtSSU=$SORTMERNADIR/rRNA_databases/mtSSU_UCLUST-95-identity.fasta

##
echo Checking

## we get two dir and two files as input
if [ $UNPAIRED == 0 ]; then
    if [ $# != 4 ]; then
	echo "This function takes two directories and two files as arguments"
	usage
    fi
else
    if [ $# != 3 ]; then
	echo "This function takes two directories and one file as argument"
	usage
    fi
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory"
    usage
fi

## 
echo Gunzipping

## unzip the files
if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing fastq.gz file"
    usage
fi
f1=`basename ${3//.gz/}`

if [ $UNPAIRED == 0 ]; then
    if [ ! -f $4 ]; then
	echo "The forth argument needs to be an existing fastq.gz file"
	usage
    fi
    f2=`basename ${4//.gz/}`
fi

## decompress them
if [ ! -f $2/$f1 ]; then
    gunzip -c $3 > $2/$f1
fi
if [ $UNPAIRED == 0 ]; then
    if [ ! -f $2/$f2 ]; then
	gunzip -c $4 > $2/$f2
    fi
fi

## interleave them
fm=`basename ${3//.f*q.gz/}`
if [ $UNPAIRED == 0 ]; then
    isVersion9=`sortmerna --version | grep "version 1.9" | wc -l`
    if [ $isVersion9 == 1 ]; then
	merge-paired-reads.sh $2/$f1 $2/$f2 $2/$fm
    else
	echo SortMeRNA must be version 1.9
	exit 1
    fi
fi

##
if [ $UNPAIRED == 0 ]; then
    echo Pre-cleaning
    rm -f $2/$f1 $2/$f2
fi

##
echo Sorting

## PE
if [ $UNPAIRED == 0 ]; then
    fo=`basename ${3//_[1,2].f*q.gz/_sortmerna}`
else
    fo=`basename ${3//.f*q.gz/_sortmerna}`
fi

## check the options
opt=
if [ $KEEP -eq 1 ]; then
    opt="--other $2/$fo"
fi 

if [ $UNPAIRED == 0 ]; then
    sortmerna -n 1 --db $mtSSU --I $2/$fm --accept $2/${fo}_mtSSU --log $1/${fo}_mtSSU -a $PROC -v --paired-in $opt
else
    sortmerna -n 1 --db $mtSSU --I $2/$f1 --accept $2/${fo}_mtSSU $1/$fo --log $1/${fo}_mtSSU -a $PROC -v $opt
fi

## deinterleave it
echo Post-Cleaning
if [ $UNPAIRED == 0 ]; then
    find $2 -name "${fo}_mtSSU*" -print0 | xargs -0 -I {} -P 6 sh -c 'unmerge-paired-reads.sh $0 $1/`basename ${0//.fastq/_1.fq}` $1/`basename ${0//.fastq/_2.fq}`' {} $1
fi

## deinterleave the rest if needed
if [ $KEEP -eq 1 ]; then
    if [ $UNPAIRED == 0 ]; then
	unmerge-paired-reads.sh $2/$fo.fastq $1/${fo}_1.fq $1/${fo}_2.fq
    fi
fi

## rm the tmp
if [ $UNPAIRED == 0 ]; then
    rm -f $2/$fm $2/$fo.fastq
else
    rm -f $2/$f1
fi

## 
echo Gzipping

## compress the output files
find $1 -name "${fo}*.fq" -print0 | xargs -0 -I {} -P 8 gzip -f {}
#printf "%s\0%s" $1/${fo}_1.fq $1/${fo}_2.fq | xargs -0 -I {} -P 2 gzip -f {}

##
echo Done
