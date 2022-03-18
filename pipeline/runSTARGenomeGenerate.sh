#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 12:00:00
#SBATCH -n 20
#SBATCH --mail-type=END,FAIL

# stop on error and undefined vars
set -eu

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# options 
BITS=18
FORMAT="gff3"
MATE=100
MEM=120000000000
THREAD=20
SAINX=14
SPARSE=1
TAG="Parent"
NOGFF=0


# USAGE
USAGETXT=\
"
	Usage: $0 [options] <singularity> <index dir> <fasta file> <gff3 file>

        Options:
            -b      the chromosome bin nbits ($BITS)
            -f      the format: gff3 or gtf ($FORMAT)
            -l      the mate length ($MATE)
            -n      no gff
            -m      maximum memory to allocate (in bytes $MEM)
            -p      number of threads to use ($THREAD)
            -s      the sparsity of the index ($SPARSE)
            -S      the SA index Nbases ($SAINX)
            -t      the tag to use in the gff3 file to link exon <-> transcript ($TAG)

       Note: if the gtf format is selected, the t option default is switched to 'transcript_id'
"
## process the options
while getopts "b:f:l:m:np:s:S:t:" opt; do
  case "$opt" in
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
if [ $# != 4 ] && [ $NOGFF -eq 0 ] ; then
    abort "This function takes a singularity container, one directory, one fasta and one gff3 file as arguments"
fi

if [ $# != 3 ] && [ $NOGFF -eq 1 ] ; then
    abort "This function takes a singularity container, one directory, one fasta file as arguments when -n is set"
fi

[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity container"

## enforce singularity
[[ -z $SINGULARITY_BINDPATH ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

[[ ! -d $2 ]] && abort "The second argument needs to be the STAR index directory"

[[ ! -f $3 ]] && abort "The third argument needs to be a fasta file"

[[ ! -f $4 ]] && [[ $NOGFF -eq 0 ]] && abort "The fourth argument needs to be a gff3 file"

## run
echo Indexing
if [ "$NOGFF" -eq "0" ]; then
  OPTIONS="--sjdbGTFfile $4 --sjdbGTFtagExonParentTranscript $TAG --sjdbOverhang $MATE"
else
  OPTIONS=
fi

singularity exec $1 STAR --runMode genomeGenerate --genomeDir $2 --genomeFastaFiles $3 \
--limitGenomeGenerateRAM $MEM --runThreadN $THREAD --genomeChrBinNbits $BITS \
--genomeSAsparseD $SPARSE --genomeSAindexNbases $SAINX $OPTIONS

## fix permission
chmod -R g+w $2

##
echo Done


