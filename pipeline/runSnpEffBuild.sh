#! /bin/bash
#SBATCH -p core -n 1
#SBATCH --mem=64GB
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
#SBATCH -J snpEff-build

# stop on error
set -ex

# usage
export USAGETXT=\
"
Usage: $0 <CONFIG-FILE> <GENOME-VERSION> <GENOME-FASTA> <GFF3> <PROTEIN-FASTA> <CDS-FASTA>
  CONFIG-FILE is the snpEff config file - the data.dir needs to be adapted and
    an entry with the genome version needs to be present: e.g.
    'sdl16.genome : Saccharomycodes_ludwigii V16'
  GENOME-VERSION is the version referenced in the snpEff config file, e.g. sdl16
  GENOME-FASTA is the genome fasta file
  GFF3 is the gene gff3 file
  PROTEIN-FASTA is the protein fasta file
  CDS-FASTA is the CDS fasta file
  
  The output will be written in the config file directory
"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# checks
if [ $# -ne 6 ]; then
  abort "This script expects 6 arguments."
fi

if [ ! -f $1 ]; then
  abort "The config file does not exist"
fi

if [ ! -f $3 ]; then
  abort "The genome fasta file does not exist"
fi

if [ ! -f $4 ]; then
  abort "The gene gff3 file does not exist"
fi

if [ ! -f $5 ]; then
  abort "The protein fasta file does not exist"
fi

if [ ! -f $6 ]; then
  abort "The CDS fasta file does not exist"
fi

# prep
out=$(dirname $1)
mkdir -p $out/$2
cd $out/$2
ln -sf $4 genes.gff
ln -sf $5 protein.fa
ln -sf $6 cds.fa
cd $out
mkdir -p $out/genomes
cd $out/genomes
ln -sf $3 $2.fa
cd $out

# run
java -jar $CLASSPATH/snpEff.jar build -c $1 $2
