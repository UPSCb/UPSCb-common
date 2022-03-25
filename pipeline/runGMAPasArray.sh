#!/bin/bash -l
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
## time too for large files
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
<<<<<<< HEAD
#SBACTH -e gmap.SLURM_ARRAY_TASK_ID.err
#SBACTH -o gmap.SLURM_ARRAY_TASK_ID.out
=======
#SBATCH -e gmap.%a.err
#SBATCH -o gmap.%a.out
>>>>>>> kogia
## mail-user and A have to be set in the submit script

## stop on error
set -e

## be verbose and extend the commands
set -x

## load the modules
#module load bioinfo-tools
#module load gmap-gsnap

## check the options if any
PROC=20
INTRONLENGTH=11000
FMT="gff3_gene"
PATHS=5
OPTIONS="-x 100 --min-identity 0.7"
CMD="gmap"

## ARRAY number
ARRAYID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

## usage
usage(){
echo >&2 \
"
	Usage: runGMAP.sh [options] <fasta file> <index dir> <index name> <out dir>
	
	Options:
                -c use cross species preset 
                -f the output format (default to gff3_gene (or $FMT))
                -i max intron length (-K and -w gmap parameters) (default $INTRONLENGTH)
                -I set min identity (default 0.7)
		-l use gmapl instead of gmap
                -n the number of paths to show (default $PATHS)
                -p number of threads to use (default $PROC)

        Note:
             You need to set the UPSCb env. variable to your UPSCb git checkout directory.
"
	exit 1
}

## get the options
while getopts m:f:i:ln:p:cI: option
do
        case "$option" in
	    f) FMT=$OPTARG;;
	    i) INTRONLENGTH=$OPTARG;;
	l) CMD="gmapl";;
        n) PATHS=$OPTARG;;
	    p) PROC=$OPTARG;;
	    c) OPTIONS="$OPTIONS --cross-species";;
		I) OPTIONS="-x 100 --min-identity "$OPTARG;;
        \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# != 4 ]; then
    echo "This function takes one file, two directories and one name as arguments"
    usage
fi

if [ ! -f $1.$ARRAYID ]; then
    echo "The first argument needs to be a fasta file"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be the GMAP index directory"
    usage
fi

if [ ! -d $2/$3 ]; then
    echo "The third argument needs to be the name of the GMAP index"
    usage
fi

if [ ! -d $4 ]; then
    echo "The forth argument needs to be the output directory"
    usage
fi

if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## run GMAP
echo Aligning

## Output format
sfx=
case $FMT in
    sampe) sfx=".sam";;
    samse) sfx=".sam";;
    gff3_gene) sfx=".gff3";;
    *) sfx=".$FMT";;
esac

## check if the file is zipped
cd $4
#case "${1: -3}" in
#    ".gz")
#	fnam=$3-`basename ${1//.fa*.gz/$sfx}`  
#	zcat $1 | $CMD -D $2 -d $3 --intronlength $INTRONLENGTH -B 5 -w $INTRONLENGTH -t $PROC -O -f $FMT -Y -n $PATHS $OPTIONS --split-output ${fnam}
#	;;
#    *)
	fnam=$3-$(basename $1.$ARRAYID)
	$CMD -D $2 -d $3 --intronlength $INTRONLENGTH -B 5 -w $INTRONLENGTH -t $PROC -O -f $FMT \
	-Y -n $PATHS $OPTIONS --split-output ${fnam} --use-shared-memory=1 $1.$ARRAYID
#	;;
#esac

##
echo Done

