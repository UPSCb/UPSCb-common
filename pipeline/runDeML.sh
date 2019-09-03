#!/bin/bash -l

#SBATCH --mem=16G
#SBATCH -t 2-00:00:00

set -ex

# default env
nULIMIT=2048
maxSample=`expr $(expr $nULIMIT / 8) - 1`
revComp="-r"

## usage
usage(){
echo >&2 \
"
	Usage: $0 [Options] <forward fastq> <reverse fastq> <forward index> <reverse index> <map file> <outfile prefix>[--] [additional deML arguments]

	Options:
                -h display this help and exit
                -f do not reverse complement I1 (forward index)
                -n the number of supported file handles (2048)
                
	Notes:
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - deML arguments.
		
		deML opens 8+ file handles per sample. If this overrides the number of file handles allowed on the
		system (2048 is the default, see the -n optoin), the map file is going to be split, consequence of 
		which, the 'unknown' reads will be removed.
"
	exit 1
}

# module load 
if [ ! -z $MODULE_VERSION ]; then
  module load stripBadStuff R bioinfo-tools deML
fi

# check the command line for the tell-tale flags
while getopts fhn: option
do
  case "$option" in
      f) revComp="";;
      h) usage;;
      n) nULIMIT=$OPTARG;;
      ?) usage;;
  esac
done
shift `expr $OPTIND - 1`

# check env
if [ -z $UPSCb ]; then
  echo "The UPSCb env. var. needs to be set to your UPSCb git checkout"
  usage
fi

# check args
if [ $# -lt 6 ]; then
  usage
fi

if [ ! -f $1 ]; then
  echo "The first argument needs to be an existing forward fastq (R1) file"
  usage
fi
R1=$1
shift

if [ ! -f $1 ]; then
  echo "The second argument needs to be an existing reverse fastq (R2) file"
  usage
fi
R2=$1
shift

if [ ! -f $1 ]; then
  echo "The third argument needs to be an existing forward index (I1) file"
  usage
fi
I1=$1
shift

if [ ! -f $1 ]; then
  echo "The fourth argument needs to be an existing reverse index (I2) file"
  usage
fi
I2=$1
shift

if [ ! -f $1 ]; then
  echo "The fifth argument needs to be an existing map file"
  usage
fi
MAP=$1
shift

if [ ! -d $(dirname $1) ]; then
  echo "The sixth argument needs to have a path pointing to an existing directory"
  usage
fi
OUT=$1
OUTDIR=$(dirname $OUT)
shift
    
# map, check the size
numSample=$(expr $(wc -l $MAP | cut -d ' ' -f 1) - 1)
    
# max up the ulimit no matter what
ulimit -n $nULIMIT

# drop the -- if provided
if [ $# != 0 ]; then
	## drop the --
	shift
fi

#strip maps of potentially incorporated spaces
tmpf=$(tempfile)
stripBadStuff -f $MAP > $tmpf
mv $tmpf $MAP

# no revcomp and small map
if [ $numSample -lt $maxSample ]; then
  deML -i $MAP -f $R1 -r $R2 -if1 $I1 -if2 $I2 -o $OUT $@  
else

  # split (and revcomp if needed)
  Rscript $UPSCb/src/R/createMapFileForDeML.R $revComp -f $MAP -o $OUTDIR -c 255
  
  # run the different maps
  for f in `find $OUTDIR -name "*_deML.txt"`; do
    deML -i $f -f $R1 -r $R2 -if1 $I1 -if2 $I2 -o $OUT $@
  done
fi  
  
# clean up
rm -f ${OUT}_unknown_[r,i][1,2].f*q.gz

# combine
for f in $(find $OUTDIR -name "*_r[1,2].fq.gz"); do
  mv $f ${f//.fq/.pass.fq}
  cat ${f//.fq/.pass.fq} ${f//.fq/.fail.fq} > $f
done
