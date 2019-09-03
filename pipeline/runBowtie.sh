#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 0-02:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

## exec
BOWTIE=

### tool sanity
if [ ! -z $SLURM_SUBMIT_DIR ]; then
    module load bioinfo-tools
    module load samtools/0.1.19
    module load bowtie/1.0.0
    BOWTIE=`which bowtie`
else
	BOWTIE=`which bowtie`
	if [ $? != 0 ]; then
		echo "please install bowtie before running this script or add it to your PATH"
		exit 1
	fi

	if [ ! -f $BOWTIE -a ! -x $BOWTIE ]; then
		echo "your bowtie does not appear to be an executable file"
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
SINGLE=0
PROC=8
## usage
usage(){
echo >&2 \
"
	Usage: runBowtie.sh [option] <genome dir> <fwd file> <rv file> <out dir> <tmp dir> [--] [additional bowtie arguments]
	
	Options:
                -e bowtie executable
                -p number of threads to be used (default: 8)
		-s if there is no reverse file 
	
	Notes:
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - bowtie arguments.

                bowtie is run with the following parameter by default:
                 -v 2 (end to end alignment with 2MMs max)
                 --best (best alignment is looked for)
                 --strata (but only in the best strata)
                 -q (input is fastq)
                 -y (try hard)
                 -m 1 (only uniquely aligning reads are reported)
                 -S (produces SAM)

"
	exit 1
}

## get the options
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

## check the arguments
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

if [ ! -f $1.1.ebwt ]; then
        echo "The genome index: $1.1.ebwt does not exists"
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
fi

## output name
outName=$output/`basename ${in1//$FIND/}`

## start star
if [ $SINGLE == 1 ]; then
    gunzip -c $in1 | $BOWTIE $@ --best --strata -m 1 -y -v 2 -q -S -p $PROC $genome - | samtools view -bS - | samtools sort - $outName
else
    if [ ! -f $tmpDir/$in1 ]; then
	gunzip -c $in1 > $tmpDir/${in1//.gz/}
    fi
    if [ ! -f $tmpDir/$in2 ]; then
	gunzip -c $in2 > $tmpDir/${in2//.gz/}
    fi

    $BOWTIE $@ --best --strata -m 1 -y -v 2 -q -S -p $PROC $genome -1 ${in1//.gz/} -2 ${in2//.gz/} | samtools view -bS - | samtools sort - $outName

    rm $tmpDir/${in1//.gz/} $tmpDir/${in2//.gz/}
 
fi

## convert sam to bam
samtools index $outName.bam
