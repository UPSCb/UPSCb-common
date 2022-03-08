#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL

# stop on error and undefined vars
set -eu

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## vars
INTRONMAX=70000
# Spruce_intronmax=70000
# Aspen_intronmax=11000
GFF=
SINGLE=0
PROC=20
FORMAT="gtf"
LIMIT=10000000000
QUANT=0
WIGGLE=0
NoGZ="--readFilesCommand zcat"

## additional options for STAR
OPTIONS="--outSAMstrandField intronMotif --outSAMmapqUnique 254 --outFilterMultimapNmax 100 \
--outReadsUnmapped Fastx --chimSegmentMin 1 --outSAMtype BAM SortedByCoordinate"

## usage
USAGETXT=\
"Usage:
    Paired end: $0 [option] <star singularity container> <samtools singularity container> <out dir> <genome dir> <genome fasta> <fwd file> <rv file> [--] [additional STAR arguments]
    Single end: $0 [option] -s <star singularity container> <samtools singularity container> <out dir> <genome dir> <genome fasta> <fastq file> [--] [additional STAR arguments]

	Options:
		-f the gtf/gff3 file format (default gtf)
		-g the path to a gtf/gff3 file
		-t quantify the transcriptome
		-l the BAM sorting memory limit ($LIMIT)
		-m the max intron length ($INTRONMAX)
		-p number of threads to be used (default: 16)
		-q set for Illumina +64 Phred score
		-s if there is no reverse
		-n no default option
		-w create wiggle files
		-z input is not compressed (default is compressed)

	Notes:
		The number of arguments is only 3 when -s is set.
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - STAR arguments.
		When the format is gff3, the exon-transcript relationship assumes a 'Parent' keylink.
"

## get the options
while getopts f:g:l:m:np:qstwz option
do
  case "$option" in
	    f) FORMAT=$OPTARG;;
	    g) GFF=$OPTARG;;
	    l) LIMIT=$OPTARG;;
	    m) INTRONMAX=$OPTARG;;
	    n) OPTIONS="";;
	    p) PROC=$OPTARG;;
	    q) OPTIONS="$OPTIONS --outQSconversionAdd -31";;
	    s) SINGLE=1;;
	    t) OPTIONS="$OPTIONS --quantMode TranscriptomeSAM"
	       QUANT=1;;
	    w) OPTION="$OPTIONS --outWigType bedGraph"
		     WIGGLE=1;;
		  z) NoGZ="";;
	    \?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

## update the options
OPTIONS="$OPTIONS $NoGZ"

## dirty if loop to accomodate for v2.3.*
if [ "$OPTIONS" != "" ]; then
  OPTIONS="$OPTIONS --limitBAMsortRAM $LIMIT"
fi

## check the arguments
echo "Parsing the arguments"
ARGS=7
if [ $SINGLE == 1 ]; then
    let "ARGS = $ARGS - 1"
    FIND=".f*.gz"
else
    FIND="_[1,2].f*q.gz"
fi

## check the number of args
[[ $# -lt $ARGS ]] && abort "This script needs 6 or 7 arguments for SE or PE data, respectively."

## get the containers
star=$1
shift
[[ ! -f $star ]] && "The first argument needs to be an existing STAR singularity container"

samtools=$1
shift
[[ ! -f $samtools ]] && "The first argument needs to be an existing samtools singularity container"

## get the out dir
outdir=$1
shift
[[ ! -d $outdir ]] && abort "The output directory: $outdir does not exist"

## get the index
genome=$1
shift
[[ ! -d $genome ]] && abort "The genome directory: $genome does not exist"

## get the genome fasta file
gfasta=$1
shift
[[ ! -f $gfasta ]] && abort "The genome fasta file: $gfasta does not exist"

## Check if the first file exists
fwd=$1
shift
[[ ! -f $fwd ]] && "The forward fastq file: $fwd does not exist"

## Check if the second file exists
[[ $SINGLE == 0 ]] && [[ ! -f $1 ]] && abort "The reverse fastq file: $1 does not exist"
[[ $SINGLE == 0 ]] && rev=$1 && shift

## if gff is set check if it exists
[[ ! -z $GFF ]] && [[ ! -f $GFF ]] && abort "The gene model gtf/gff3 file: $GFF does not exists"
[[ ! -z $GFF ]] && OPTIONS="--sjdbGTFfile $GFF $OPTIONS"

## if format is set
case $FORMAT in
    gff3)
    OPTIONS=" $OPTIONS --sjdbGTFtagExonParentTranscript Parent"
    ;;
    gff)
    OPTIONS=" $OPTIONS --sjdbGTFtagExonParentTranscript Parent"
    ;;
    gtf);;
    #nothing to do
    *)
	abort "There are only 2 supported format, gtf or gff3"
esac

## do we have more arguments? drop the --
[[ $# != 0 ]] && shift

## output prefix
bnam=$(basename ${fwd//$FIND/})
fnam=$outdir/$bnam

## start STAR
echo "Aligning"
if [ $SINGLE == 1 ]; then
    singularity exec $star STAR --genomeDir $genome --readFilesIn $fwd --runThreadN $PROC \
    --alignIntronMax $INTRONMAX --outFileNamePrefix $fnam $OPTIONS $@
else
    singularity exec $star STAR --genomeDir $genome --readFilesIn $fwd $rev --runThreadN $PROC \
    --alignIntronMax $INTRONMAX --outFileNamePrefix $fnam $OPTIONS $@
fi

## save the log
echo "Logging"
mkdir -p ${fnam}_logs
mv ${fnam}Log.* ${fnam}_logs

## save the junctions
mkdir -p ${fnam}_junctions
mv ${fnam}SJ* ${fnam}_junctions
mv ${fnam}Chimeric.out.junction ${fnam}_junctions

## save the wig
if [ $WIGGLE -eq 1 ]; then
	echo "Wiggling"
	mkdir -p ${fnam}_bedgraphs
	mv ${fnam}Signal.*.bg ${fnam}_bedgraphs
fi

## rename the output
echo "Renaming"
mv ${fnam}Aligned.sortedByCoord.out.bam ${fnam}_STAR.bam
if [ $SINGLE == 0 ]; then
    mv ${fnam}Unmapped.out.mate1 ${fnam}_Unmapped_1.fq
    mv ${fnam}Unmapped.out.mate2 ${fnam}_Unmapped_2.fq
else
    mv ${fnam}Unmapped.out.mate1 ${fnam}_Unmapped.fq
fi

## compress files (we would only need 2 CPUS, but what if PROC is set to 1)
find $outdir -name "${bnam}_Unmapped*.fq" -print0 | xargs -P $PROC -0 -I {} gzip -f {}

## sort the transcriptome bam and rename
if [ $QUANT == 1 ]; then
  mv ${fnam}Aligned.toTranscriptome.out.bam ${fnam}_STAR_Transcriptome.bam
  singularity exec $samtools samtools sort -@ 16 -n ${fnam}_STAR_Transcriptome.bam -o ${fnam}_STAR_Transcriptome.sorted.bam
  rm ${fnam}_STAR_Transcriptome.bam
  mv ${fnam}_STAR_Transcriptome.sorted.bam ${fnam}_STAR_Transcriptome.bam
fi

## convert the chimeric sam to cram
singularity exec $samtools samtools view -CT $gfasta ${fnam}Chimeric.out.sam | \
singularity exec $samtools samtools sort -@ 16 - -o ${fnam}_STAR_Chimeric.cram

## convert the output BAM in CRAM
singularity exec $samtools samtools view -CT $gfasta -o ${fnam}_STAR.cram ${fnam}_STAR.bam

## index the CRAMs
echo "Indexing"
printf "%s\0%s" ${fnam}_STAR.cram ${fnam}_STAR_Chimeric.cram | \
xargs -P $PROC -0 -I {} singularity exec $samtools samtools index {}

## cleanup
echo "Cleaning"
rm  ${fnam}Chimeric.out.sam ${fnam}_STAR.bam
rm -rf ${fnam}_STARtmp/
