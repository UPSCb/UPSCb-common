#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 0-02:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -eux

## helper function
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## Variables
BOWTIE=
SINGLE=0
PROC=20
OPT="-k 50 --non-deterministic"

## Usage
USAGETXT=\
"
	Usage: runBowtie2.sh [option] <genome dir> <fwd file> <rv file> <out dir> <tmp dir> [--] [additional bowtie arguments]

	Options:
                -e bowtie executable
                -p number of threads to be used (default: $PROC)
		        -s if there is no reverse file

	Notes:
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - bowtie arguments.

        bowtie2 is run with the following parameter by default:
                 -k 50 --non-deterministic

        Setting the [additional bowtie arguments] will override these defaults, set them
        again if you want to extend them

        The number of processes is set using the -p option of the script, do not set it as
        part of the [additional bowtie arguments]
"

## Tool sanity
if [ ! -z $SLURM_SUBMIT_DIR ]; then
    module load bioinfo-tools
    module load samtools
    module load bowtie2
    BOWTIE=`which bowtie2`
else
	BOWTIE=`which bowtie2`
	if [ $? != 0 ]; then
		abort "please install bowtie before running this script or add it to your PATH"
	fi

	if [ ! -f $BOWTIE -a ! -x $BOWTIE ]; then
		abort "your bowtie does not appear to be an executable file"
	fi

	samtools=`which samtools`
	if [ $? != 0 ]; then
		abort "please install samtools before running this script or add it to your PATH"
	fi

	if [ ! -f $samtools -a ! -x $samtools ]; then
		abort "your samtools does not appear to be an executable file"
	fi
fi

## Options
while getopts e:sp: option
do
        case "$option" in
	    e) BOWTIE=$OPTARG;;
	    p) PROC=$OPTARG;;
	    s) SINGLE=1;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## Arguments
if [ ! -z $BOWTIE -a ! -f $BOWTIE -a ! -x $BOWTIE ]; then
    echo "your bowtie does not appear to be an executable file"
    exit 1
fi

ARGS=5
if [ $SINGLE == 1 ]; then
    let "ARGS = $ARGS - 1"
    FIND=".f*q.gz"
else
    FIND="_[1,2].f*q.gz"
fi

if [ $# -lt $ARGS ]; then
    echo "This script needs 4 arguments for SE data; 5 for PE."
    usage
fi

if [ ! -f $1.1.bt2 ]; then
        echo "The genome index: $1.1.bt2 does not exists"
        usage
else
        genome=$1
        shift
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
        echo "The output directory: $1  does not exist"
        usage
else
        output=$1
        shift
fi

if [ ! -d $1 ]; then
        echo "The tmp directory: $1  does not exist"
        usage
else
        tmpDir=$1
        shift
fi

## do we have more arguments
if [ $# != 0 ]; then
	## drop the --
	shift
	OPT=
fi

## output name
outName=$output/`basename ${in1//$FIND/.bam}`

## start star
if [ $SINGLE == 1 ]; then
    gunzip -c $in1 | $BOWTIE $@ $OPT -p $PROC -x $genome -U - | samtools view -bS - | samtools sort -o $outName -
else
    if [ ! -f $tmpDir/$in1 ]; then
	    gunzip -c $in1 > $tmpDir/${in1//.gz/}
    fi
    if [ ! -f $tmpDir/$in2 ]; then
	    gunzip -c $in2 > $tmpDir/${in2//.gz/}
    fi

    $BOWTIE $@ $OPT -p $PROC $genome -1 ${in1//.gz/} -2 ${in2//.gz/} | samtools view -bS - | samtools sort - $outName

    rm $tmpDir/${in1//.gz/} $tmpDir/${in2//.gz/}
 
fi

## index
samtools index $outName

## Done
