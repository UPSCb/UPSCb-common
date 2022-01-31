#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 12:00:00
#SBATCH -n 20
#SBATCH --mail-type=ALL

## stop on error and be verbose
set -ex

## load the modules
#module load bioinfo-tools
#module load star

## usage
usage(){
echo >&2 \
"
	Usage: $0 <index dir> <fasta file> <gff3 file>

        Options:
            -b      the chromosome bin nbits (18)
            -f      the format: gff3 or gtf (gff3)
            -l      the mate length (100)
            -n      no gff
            -m      maximum memory to allocate (120GB, in bytes 120000000000)
            -p      number of threads to use (20)
            -s      the sparsity of the index (1)
            -S      the SA index Nbases (14)
            -t      the tag to use in the gff3 file to link exon <-> transcript (Parent)

       Note: if the gtf format is selected, the t option default is switched to 'transcript_id'
"
	exit 1
}

## process the options
BITS=18
FORMAT="gff3"
MATE=100
MEM=120000000000
THREAD=20
SAINX=14
SPARSE=1
TAG="Parent"
NOGFF=0

while getopts "b:f:l:m:np:sS:t:" opt; do
    case $opt in
	b) BITS=$OPTARG;;
	f) FORMAT=$OPTARG;;
	l) MATE=$OPTARG;;
	n) NOGFF=1;;
	m) MEM=$OPTARG;;
        p) THREAD=$OPTARG;;
	s) SPARSE=$OPTARG;;
	S) SAINX=$OPTARG;;
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
if [ $# != 3 ] && [ $NOGFF -eq 0 ] ; then
    echo "This function takes one directory, one fasta and one gff3 file as arguments"
    usage
fi

if [ $# != 2 ] && [ $NOGFF -eq 1 ] ; then
    echo "This function takes one directory, one fasta file as arguments when -n is set"
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

if [ ! -f $3 ] && [ $NOGFF -eq 0 ] ; then
    echo "The third argument needs to be a gff3 file"
    usage
fi

## run
echo Indexing
if [ "$NOGFF" -eq "0" ]; then
  OPTIONS="--sjdbGTFfile $3 --sjdbGTFtagExonParentTranscript $TAG --sjdbOverhang $MATE"
else
  OPTIONS=
fi

STAR --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 \
--limitGenomeGenerateRAM $MEM --runThreadN $THREAD --genomeChrBinNbits $BITS \
--genomeSAsparseD $SPARSE --genomeSAindexNbases $SAINX $OPTIONS

## fix permission
chmod -R g+w $1

##
echo Done


