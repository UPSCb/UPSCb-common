#!/bin/bash
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL

set -ex

#module load bioinfo-tools GATK

# vars
Threads=8
JavaMem=40G

## usage
usage(){
echo >&2 \
"
  Usage: $0 [options] <in.bam> <ref.fasta> <output directory> [--] [hc_args ...]
  
  Options: -m: the memory argument for java, default to 40G
           -t: the number of threads to use, default to 8
  Notes: -- is a special argument that stop the command line scanning for the script  options.      
"
exit 1
}

## get the options
while getopts t:m: option
do
  case "$option" in
	    t) Threads=$OPTARG;;
	    m) JavaMem=$OPTARG;;
		\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`


#  check
if [ ! -d "$GATK_HOME" ]; then
    echo >&2 "Could not find GATK"
    usage
fi

#setup
if [ -d "$SNIC_NOBACKUP" ]; then
    tmp=$SNIC_NOBACKUP
else
    tmp=/mnt/picea/tmp
fi

if [ $# -lt 3 ]; then
    usage
fi

if [ ! -d $tmp ]; then
    echo "tmp is not a directory" 1>&2
    usage
fi

# get multiple files
if [ $(echo $1 | grep -c ",") -eq 1 ]; then 
  echo "Checking the file list. If this fails, file(s) may be missing."
  echo $1 | xargs -d, -P $Threads -I {} bash -c 'if [ ! -f $0 ]; then exit 1; fi' {}
  inbam=$(echo $1 | sed "s:,: -I :")
  bname=$(basename $(echo $inbam | sed "s: .*::"))
else
  if [ ! -f $1 ]; then
      echo "Could not find BAM file '$1'" 1>&2
      usage
  fi
  inbam=$1
  bname=`basename $inbam`
fi

# ref
if [ ! -f $2 ]; then
    echo "Could not find reference '$1'" 1>&2
    usage
fi
ref=$2

# outdir
if [ ! -d $3 ]; then
    echo "No such output directory" 1>&2
    usage
fi
outdir=$3

# get additional args
shift 3
if [ $# != 0 ]; then
	## drop the --
	shift
fi

# file names
sname="${bname/_[st]*[st]*_STAR*.bam/}"
outfile="$outdir/${sname}.ug.raw.snps.indels.vcf"

# Run the UnifiedGenotyper
java -Xmx${JavaMem} -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -I $inbam \
    -R $ref \
    -o $outdir/$sname.ug.raw.snps.indels.vcf \
    -nct $Threads \
    -T UnifiedGenotyper \
    $@

