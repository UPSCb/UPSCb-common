#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 12:00:00
#SBATCH -n 16
#SBATCH --mem=100G
#SBATCH --mail-type=ALL

## stop on error and be verbose
set -ex

## load the modules
module load bioinfo-tools
module load star

## usage
usage(){
echo >&2 \
"
	Usage: $0 <index dir> <fasta file> <gff3 file>

        Options:
            -b      the chromosome bin nbits (18)
            -f      the format: gff3 or gtf (gff3)
            -l      the mate length (100)
            -m      maximum memory to allocate (100GB, in bytes 100000000000)
            -p      number of threads to use (16)
            -s      the sparsity of the index (1)
            -t      the tag to use in the gff3 file to link exon <-> transcript (Parent)

       Note: if the gtf format is selected, the t option default is switched to 'transcript_id'
"
	exit 1
}

## process the options
BITS=18
FORMAT="gff3"
MATE=100
MEM=100000000000
THREAD=16
SPARSE=1
TAG="Parent"
while getopts "b:f:l:m:p:s:t:" opt; do
    case $opt in
	b) BITS=$OPTARG;;
	f) FORMAT=$OPTARG;;
	l) MATE=$OPTARG;;
	m) MEM=$OPTARG;;
        p) THREAD=$OPTARG;;
	s) SPARSE=$OPTARG;;
	t) TAG=$OPTARG;;
        \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

## change the option if FORMAT is gtf and the default is set
if [ $FORMAT == "gtf" ] && [ $TAG == "Parent" ]; then
    TAG="transcript_id"
fi

## we get one out dir, one fasta and one gff3 file as input
if [ $# != 3 ]; then
    echo "This function takes one directory, one fasta and one gff3 file as arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the STAR index directory"
    usage
fi

if [ ! -f $2 ]; then
    echo "The second argument needs to be a fasta file"
    usage
fi

if [ ! -f $3 ]; then
    echo "The third argument needs to be a gff3 file"
    usage
fi

## run GMAP
echo Indexing

## run
STAR --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --sjdbGTFfile $3 --limitGenomeGenerateRAM $MEM --runThreadN $THREAD --genomeChrBinNbits $BITS --genomeSAsparseD $SPARSE --sjdbGTFtagExonParentTranscript $TAG --sjdbOverhang $MATE

## fix permission
chmod -R g+w $1

##
echo Done


