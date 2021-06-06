#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 0-12:00:00
#SBATCH --mail-type=ALL

# -p node is needed to accept the -C memory configuration

## stop on error and be verbose in the output
set -e -x

## load the modules
#module load bioinfo-tools star samtools

## vars
INTRONMAX=70000
# Spruce_intronmax=1000000
# Aspen_intronmax=11000
GFF=
SINGLE=0
PROC=20
FORMAT="gtf"
LIMIT=10000000000
QUANT=0
UNSORTED=0
BEDGRAPH=0
CHIMERIC=0
CRAM=1

## additional options for STAR
OPTIONS="--readFilesCommand zcat --outSAMmapqUnique 254 --outFilterMultimapNmax 100"
#OPTIONS="$OPTIONS --twopassMode Basic"

## usage
usage(){
echo >&2 \
"Usage:
    Paired end: $0 [option] <out dir> <genome dir> <genome fasta> <fwd file> <rv file> [--] [additional STAR arguments]
    Single end: $0 [option] -s <out dir> <genome dir> <genome fasta> <fastq file> [--] [additional STAR arguments]

	Options:
	  -b produce bedgraphs
	  -c produce chimeric files
    -f the gtf/gff3 file format (default gtf)
    -g the path to a gtf/gff3 file
    -h print the usage
    -i add intronMotif as a SAM attribute
		-l the BAM sorting memory limit ($LIMIT)
		-m the max intron length ($INTRONMAX)
		-n no default option
		-o simple bam unsorted output
    -p number of threads to be used (default: $PROC)
		-q set for Illumina +64 Phred score
		-s if there are no reverse reads (single-end mode)
		-t quantify the transcriptome
		-z do not bam to cram

	Notes:
		The number of arguments is only 3 when -s is set.
		-- is a special argument that stop the command line scanning for the script  options.
		It is necessary if you want to precised additional - non-default - STAR arguments.
		When the format is gff3, the exon-transcript relationship assumes a 'Parent' keylink.
"
	exit 1
}

## get the options
while getopts bcf:g:hil:m:nop:qstz option
do
  case "$option" in
      b) OPTIONS="$OPTIONS --outWigType bedGraph"
         BEDGRAPH=1;;
      c) OPTIONS="$OPTIONS --chimSegmentMin 1"
         CHIMERIC=1;;
	    f) FORMAT=$OPTARG;;
	    g) GFF=$OPTARG;;
	    h) usage;;
	    i) OPTIONS="$OPTIONS --outSAMstrandField intronMotif";;
	    l) LIMIT=$OPTARG;;
      m) INTRONMAX=$OPTARG;;
      o) UNSORTED=1;;
      n) OPTIONS="";;
	    p) PROC=$OPTARG;;
	    q) OPTIONS="$OPTIONS --outQSconversionAdd -31";;
	    s) SINGLE=1;;
	    t) OPTIONS="$OPTIONS --quantMode TranscriptomeSAM"
	       QUANT=1;;
	    z) CRAM=0;;
	    \?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

## update the options
## dirty if loop to accomodate for v2.3.*
if [ "$OPTIONS" != "" ]; then
  OPTIONS="$OPTIONS --limitBAMsortRAM $LIMIT"
fi
OPTIONS="$OPTIONS --alignIntronMax $INTRONMAX"
OPTIONS="$OPTIONS --runThreadN $PROC"

## check the arguments
echo "Parsing the arguments"
ARGS=5
if [ $SINGLE == 1 ]; then
    let "ARGS = $ARGS - 1"
    FIND=".f*.gz"
else
    FIND="_[1,2].f*q.gz"
fi

## checkthe number of args
if [ $# -lt $ARGS ]; then
    echo "This script needs 4 or 5 arguments for SE or PE data, respectively."
    usage
fi

## get the out dir
outdir=$1
shift

## check the genome dir
if [ ! -d $1 ]; then
        echo "The genome directory: $1 does not exist"
        usage
else
    genome=$1
    shift
fi

## get the genome fasta file
if [ ! -f $1 ]; then
  echo "The genome fasta file: $1 does not exist"
	usage
else
  gfasta=$1
  shift
fi

## Check if the first file exists
if [ ! -f $1 ]; then
	echo "The forward fastq file: $1 does not exist"
	usage
else
    fwd=$1
    shift
fi

## Check if the second file exists
if [ $SINGLE == 0 ]; then
    if [ ! -f $1 ]; then
        echo "The reverse fastq file: $1 does not exist"
        usage
    else
	rev=$1
	shift
    fi
fi

## if gff is set check if it exists
if [ ! -z $GFF ] && [ ! -f $GFF ] ; then
    echo "The gene model gtf/gff3 file: $GFF does not exists"
    usage
else
  if [ ! -z $GFF ]; then
    OPTIONS="--sjdbGTFfile $GFF $OPTIONS"
  fi
fi

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
	echo "There are only 2 supported format, gtf or gff3"
	usage;;
esac

## output
if [ $UNSORTED -eq 1 ]; then
  OPTIONS="$OPTIONS --outSAMtype BAM Unsorted --outSAMattributes None"
else
  OPTIONS="$OPTIONS --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMattributes All"
fi


## do we have more arguments
if [ $# != 0 ]; then
	## drop the --
	shift
fi

## create the output dir
echo "Processing"
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

## output prefix
bnam=`basename ${fwd//$FIND/}`
fnam=$outdir/$bnam

## start STAR
echo "Aligning"
if [ $SINGLE == 1 ]; then
    STAR --genomeDir $genome --readFilesIn $fwd --outFileNamePrefix $fnam $OPTIONS $@
else
    STAR --genomeDir $genome --readFilesIn $fwd $rev --outFileNamePrefix $fnam $OPTIONS $@
fi

## save the log
echo "Logging"
mkdir -p ${fnam}_logs
mv ${fnam}Log.* ${fnam}_logs

## save the junctions
mkdir -p ${fnam}_junctions
mv ${fnam}SJ* ${fnam}_junctions
[[ $CHIMERIC -eq 1 ]] && mv ${fnam}Chimeric.out.junction ${fnam}_junctions

## save the wig
if [ $BEDGRAPH -eq 1 ]; then
  echo "Wiggling"
  mkdir -p ${fnam}_bedgraphs
  mv ${fnam}Signal.*.bg ${fnam}_bedgraphs
fi

## rename the output
echo "Renaming"
if [ $UNSORTED -eq 0 ]; then
  mv ${fnam}Aligned.sortedByCoord.out.bam ${fnam}_STAR.bam
  if [ $SINGLE == 0 ]; then
    mv ${fnam}Unmapped.out.mate1 ${fnam}_Unmapped_1.fq
    mv ${fnam}Unmapped.out.mate2 ${fnam}_Unmapped_2.fq
  else
    mv ${fnam}Unmapped.out.mate1 ${fnam}_Unmapped.fq
  fi
  ## compress files (we would only need 2 CPUS, but what if PROC is set to 1)
  find $outdir -name "${bnam}_Unmapped*.fq" -print0 | xargs -P $PROC -0 -I {} gzip -f {}
else
  mv ${fnam}Aligned.out.bam ${fnam}_STAR.bam
fi

## sort the transcriptome bam and rename
if [ $QUANT == 1 ]; then
  mv ${fnam}Aligned.toTranscriptome.out.bam ${fnam}_STAR_Transcriptome.bam
  samtools sort -@ $PROC -n ${fnam}_STAR_Transcriptome.bam -o ${fnam}_STAR_Transcriptome.sorted.bam
  rm ${fnam}_STAR_Transcriptome.bam
  mv ${fnam}_STAR_Transcriptome.sorted.bam ${fnam}_STAR_Transcriptome.bam
fi

## convert the chimeric sam to cram
if [ $CRAM -eq 1 ]; then
  if [ $CHIMERIC -eq 1 ]; then
    samtools view -CT $gfasta ${fnam}Chimeric.out.sam | samtools sort -@ $PROC - -o ${fnam}_STAR_Chimeric.cram
    samtools index ${fnam}_STAR_Chimeric.cram 
  fi
  
  ## convert the output BAM in CRAM
  samtools view -CT $gfasta -o ${fnam}_STAR.cram ${fnam}_STAR.bam

  ## index the CRAMs
  echo "Indexing"
  samtools index ${fnam}_STAR.cram

  ## cleanup
  echo "Cleaning"
  rm ${fnam}_STAR.bam
else
  if [ $CHIMERIC -eq 1 ]; then 
    samtools view -b ${fnam}Chimeric.out.sam | samtools sort -@ $PROC - -o ${fnam}_STAR_Chimeric.bam
    samtools index ${fnam}_STAR_Chimeric.bam
  fi
  echo "Indexing"
  #samtools index ${fnam}_STAR.bam

  ## cleanup
  echo "Cleaning"
fi
[[ $CHIMERIC -eq 1 ]] && rm ${fnam}Chimeric.out.sam
rm -rf ${fnam}_STARtmp/
