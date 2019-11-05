#!/bin/bash
#SBATCH -p core
#SBATCH -n 8
#SBATCH --mem=40GB
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

set -ex

# module load bioinfo-tools GATK

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# vars
#Threads=8
JavaMem=40G

USAGETXT=\
"
Usage: $0 [options] <in.bam> <ref.fasta> <output directory> [--] [hc_args ...]

Options: -m: the memory argument for java, default to 40G

Notes: -- is a special argument that stop the command line scanning for the script  options.      
      
      This script is GATK v4 compatible and GATK V3 incompatible. More at https://software.broadinstitute.org/gatk/blog?id=7847
      
      GATK 4 will determine how many CPUs are available and use these.
"

## get the options
while getopts m: option
do
  case "$option" in
	    m) JavaMem=$OPTARG;;
		\?) ## unknown flag
		usage;;
  esac
done
shift `expr $OPTIND - 1`

#  check
isExec gatk

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
    abort "tmp is not a directory"
fi

if [ $(echo $1 | grep -c ",") -eq 1 ]; then 
  echo "Checking the file list. If this fails, file(s) may be missing."
  echo $1 | xargs -d, -P $Threads -I {} bash -c 'if [ ! -f $0 ]; then exit 1; fi' {}
  inbam=$(echo $1 | sed "s:,: -I :")
  bname=$(basename $(echo $inbam | sed "s: .*::"))
else
  if [ ! -f $1 ]; then
      abort "Could not find BAM file '$1'"
  fi
  inbam=$1
  bname=`basename $inbam`
fi

if [ ! -f $2 ]; then
    abort "Could not find reference '$2'"
fi
ref=$2

if [ ! -d $3 ]; then
    abort "No such output directory"
fi
outdir=$3

# get additional args
shift 3
if [ $# != 0 ]; then
	## drop the --
	shift
fi

sname="${bname/_[st]*[st]*_STAR*.bam/}"
outfile="$outdir/${sname}_raw.snps.indels.vcf"

# Run the HaplotypeCaller
gatk --java-options "-Xmx${JavaMem} -Djava.io.tmpdir=$tmp" HaplotypeCaller -I $inbam -R $ref -O $outfile  $@ 
