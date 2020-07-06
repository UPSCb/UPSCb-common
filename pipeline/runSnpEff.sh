#! /bin/bash
#SBATCH -p core -n1
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL
#SBATCH -J snpEff-build
#SBATCH --mem=16G

# stop on error
set -ex

# usage
export USAGETXT=\
"
Usage: $0 <CONFIG-FILE> <GENOME-VERSION> <VCF>
  CONFIG-FILE is the snpEff config file - the data.dir needs to be adapted and
    an entry with the genome version needs to be present: e.g.
    'sdl16.genome : Saccharomycodes_ludwigii V16'
  GENOME-VERSION is the version referenced in the snpEff config file, e.g. sdl16
  VCF is the vcf file to process

  The output will be written in the vcf input file directory. Best is to link the
  input vcf file there
"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# checks
if [ $@ -ne 3 ]; then
  abort "This script expects 3 arguments."
fi

if [ ! -f $1 ]; then
  abort "The config file does not exist"
fi

if [ ! -f $3 ]; then
  abort "The vcf file does not exist"
fi

# prep
out=$(dirname $3)
fnam=$(basename ${3/.vcf.*/})

# run
java -Xms4g -Xmx16g -jar $CLASSPATH/snpEff.jar ann -c $1 \
-s $out/${fnam}_snpEff_summary.html $2 $3 > $out/$fnam.ann.vcf

# compress and index
bgzip -f $out/$fnam.ann.vcf
tabix $out/$fnam.ann.vcf.gz -p vcf
