#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 60
#SBATCH -t 15-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -w watson

## stop on error
set -e
set -x

module load bioinfo-tools framedp/1.2.2

## defaults
CFG=$3

# fetches a --notrain argument
options=$5

ref_folder=$4

## usage
usage(){
echo >&2 \
"
        Usage: `basename $0` <fasta file> <out dir> <CFG> <ref path> options
        
        Note: if the fasta file is gzipped, it will be unzipped in the out
        directory temporarily
"
 exit 1
}

## check
if [ $# -gt 5 ] || [ $# -lt 4 ]; then
	echo "Not correct number of arguments
	"
	echo $#
    usage
fi

if [ ! -f $1 ] ; then
    echo "The first argument should be the fasta file." $1
    usage
fi

if [ ! -d $2 ] ; then
    echo "The second argument should be a folder" $2
    usage
fi

if [ ! -f $3 ] ; then
    echo "The third argument has to be a CFG file" $3
    usage
fi

if [ ! -d $4 ] ; then
    echo "The forth argument has to be a reference directory" $4
    usage
fi

if [ ! -d $2/workdir ] ; then
    mkdir -p $2
    mkdir -p $2/workdir
fi

## In case of a -notrain run Create local reference folder, copy trainingset from reference folder 
## linking is not enough as we need to overwrite some files
if [ "$ref_folder" != "" ] && [ "$options" == "--no_train" ]; then
    cp -r $ref_folder/* $2
fi

## If the input file is zipped, decompress
input=`readlink -f $1`
clean=0
if [ "${input##*.}" == "gz" ]; then
  gunzip -c $input > $2/`basename ${input//.gz/}`
  input=$2/`basename ${input//.gz/}`
  clean=1
fi

##run
$FRAMEDP/bin/FrameDP.pl --cfg $CFG --infile $input --outdir `readlink -f $2` --workingdir `readlink -f $2`/workdir $options

if [ $clean == 1 ]; then
  rm $input
fi

#if [ $options == "--notrain" ] && [ $ref_folder != "" ]; then
#    cp -r $2 $ref_folder
#fi



