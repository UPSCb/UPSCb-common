#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 0-06:00:00
#SBATCH --mail-type=ALL

#################
## Build geneModel
#################
## TODO extract that to its own script
##usage  sbatch -p devel -t 1:00:00 runSTAR.sh genome.fa
#/home/davidsu/bin/STAR --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --sjdbOverhang 99 --sjdbGTFfile $3 --runThreadN 8
#exit;

## stop on error and be verbose in the output
set -e -x


## exec
STAR=

### tool sanity
if [ ! -z $SLURM_SUBMIT_DIR ]; then
    module load bioinfo-tools
    module load samtools/0.1.19
    module load star/2.3.0e
    STAR=`which STAR`
else
	STAR=`which STAR`
	if [ $? != 0 ]; then
		echo "please install STAR before running this script or add it to your PATH"
		exit 1
	fi

	if [ ! -f $STAR -a ! -x $STAR ]; then
		echo "your STAR does not appear to be an executable file"
		exit 1
	fi

	samtools=`which samtools`
	if [ $? != 0 ]; then
		echo "please install samtools before running this script or add it to your PATH"
		exit 1
	fi

	if [ ! -f $samtools -a ! -x $samtools ]; then
		echo "your samtools does not appear to be an executable file"
		exit 1
	fi
fi

##########
# Run star
##########
INTRONMAX=11000
OUT_DIR=`pwd`
GFF=1
SINGLE=0

## usage
usage(){
echo >&2 \
"
	Usage: runSTAR.sh [option] <fwd file> <rv file> <genome dir> <gene model gff3> [--] [additional STAR arguments]
	
	Options:
                -e STAR executable
                -g if there is no gff file
		-m  max intron length
		-o  outdir
                -s if there is no reverse file 
	
	Notes:
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - STAR arguments.
"
	exit 1
}

## get the options
while getopts e:gm:o:s option
do
        case "$option" in
	    e) STAR=$OPTARG;;
	    g) GFF=0;;
            m) INTRONMAX=$OPTARG;;
            o) OUT_DIR=$OPTARG;;
	    s) SINGLE=1;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## check the arguments
if [ ! -z $STAR -a ! -f $STAR -a ! -x $STAR ]; then
    echo "your STAR does not appear to be an executable file"
    exit 1
fi

ARGS=4
if [ $SINGLE == 1 ]; then
    let "ARGS = $ARGS - 1"
    FIND=".f*q.gz"
else
    FIND="_1.f*q"
fi

if [ $GFF == 0 ]; then
    let "ARGS = $ARGS - 1"
fi

if [ $# -lt $ARGS ]; then
    echo "This script needs 2 arguments without GFF and for SE data; 3 for either and 4 for none of these two conditions."
    usage
fi

if [ ! -f $1 ]; then
	echo "The forward fastq file: $1 does not exist"
	usage	
else
	in1=$1
	shift
fi

if [ $SINGLE == 0 ]; then
    if [ ! -f $1 ]; then
        echo "The reverse fastq file: $1 does not exist"
        usage
    else
        in2=$1
        shift
    fi
fi

if [ ! -d $1 ]; then
        echo "The genome directory: $1  does not exist"
        usage
else
        genome=$1
        shift
fi

if [ $GFF == 1 ]; then
    if [ ! -f $1 ]; then
        echo "The gene model gff3 file: $1  does not exists"
        usage
    else
        gff3=$1
        shift
    fi
fi

## do we have more arguments
if [ $# != 0 ]; then
	## drop the --
	shift
fi

## output name
uz3=$OUT_DIR/`basename ${in1//$FIND/}`

## start star
if [ $SINGLE == 1 -a $GFF == 0 ]; then
    $STAR --genomeDir $genome --readFilesIn $in1 --runThreadN 8 --alignIntronMax $INTRONMAX --outSAMstrandField intronMotif --readFilesCommand zcat --outFileNamePrefix $uz3 $@
else
    if [ $SINGLE == 1 -o $GFF == 0 ]; then
	if [ $GFF == 0 ]; then
	    $STAR --genomeDir $genome --readFilesIn $in1 $in2 --runThreadN 8 --alignIntronMax $INTRONMAX --outSAMstrandField intronMotif --outFileNamePrefix $uz3 $@
	else
	    $STAR --genomeDir $genome --readFilesIn $in1 --runThreadN 8 --alignIntronMax $INTRONMAX --outSAMstrandField intronMotif --sjdbGTFfile $gff3 --readFilesCommand zcat --outFileNamePrefix $uz3 $@
	fi
    else
	$STAR --genomeDir $genome --readFilesIn $in1 $in2 --runThreadN 8 --alignIntronMax $INTRONMAX --outSAMstrandField intronMotif --sjdbGTFfile $gff3 --readFilesCommand zcat --outFileNamePrefix $uz3 $@
    fi
fi

## save the logs
mkdir -p ${uz3}_logs
mv ${uz3}Log.* ${uz3}_logs
mv ${uz3}SJ* ${uz3}_logs

## convert sam to bam
samtools view -Sb ${uz3}Aligned.out.sam | samtools sort - ${uz3}_STAR
samtools index ${uz3}_STAR.bam

## clean
rm ${uz3}Aligned.out.sam

