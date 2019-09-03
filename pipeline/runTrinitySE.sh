#!/bin/bash -l
#SBATCH -p core
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 16
## time too for large files
#SBATCH --mail-type=ALL

## stop on error, be verbose and expand the commands
set -e -x

## load necessary modules
#module load bioinfo-tools
#module load trinity

## check the options if any
MEM=40G
KMER=2
PROC=16
NORM=
BOWTIE=

## usage
usage(){
echo >&2 \
"
	Usage: runTrinity.sh [options] <out dir> <left fq(s)> <right fq(s)> <single fq(s)>
	
	Options:
                -b do not run bowtie as part of Chrysalis. The only reason to do that
                   is if the inchworm file is larger than 4Gbp. bowtie cannot index
                   files larger than 4Gbp. It is also extremely slow (samtools sort 
                   step is mono-CPU)
                -k min kmer cov (default 2 w/o digital normalization, 1 with) 
                -m mem requirement (default 40G) (rule of thumb 1.5G/1M read pairs)
                -n digitally normalize the data (coverage 50X)
		-p number of threads to use (default 16)
		-s strand specific RNAseq read orientation if paired: RF or FR, if single F or R 
		   dUTP method=RF 

        Note:
             Several fastq files can be provided as argument; separated by commas
             The 4th argument (the single fq(s)) is optional and can be left empty
"
	exit 1
}

## get the options
while getopts bk:m:p:s: option
do
        case "$option" in
	    b) BOWTIE="--no_bowtie";;
	    k) KMER=$OPTARG;;
	    m) MEM=$OPTARG;;
	    p) PROC=$OPTARG;;
	    s) SPEC="--SS_lib_type $OPTARG";;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## change default
if [ ! -z $NORM ] && [ $KMER -eq 2 ]; then
    KMER=1
fi

## we get one dir and one file as input
if [ $# -lt 2 ] ; then
    echo "This function needs 2 arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument (output dir) needs to be an existing directory"
    usage
fi

## proceed the other args
echo "Checking the fastq file(s). If this fails, file(s) may be missing."
echo $2 | xargs -d, -P $PROC -I {} bash -c 'if [ ! -f $0 ]; then exit 1; fi' {}

## nothing to do, trinity handles gz files
#echo "Preparing fwd and rev"

## add the single to the command line
#SGL=
#if [ ! -z $4 ]; then
#    echo Adding single files
#    SGL="--single `echo $4 | tr ',' ' '`"
#fi

## if there are too many files, let's do the diginorm in chunk
if [ `echo $2 | grep -o "," | wc -l` -gt 50 ] && [ ! -z $NORM ]; then
#    NORM="$NORM --normalize_by_read_set"
    NORM="--normalize_by_read_set"
fi

## run trinity
echo Assembling

singularity exec --bind /mnt:/mnt /mnt/picea/projects/singularity/trinity-2.8.3.1.simg \
Trinity  $BOWTIE --seqType fq --max_memory $MEM --single `echo $2 | tr ',' ' '` \
--output $1 --CPU $PROC $NORM $SPEC

## --min_kmer_cov $KMER 
echo Done

