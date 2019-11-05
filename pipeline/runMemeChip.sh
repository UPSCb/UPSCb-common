#! /bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL
set -eux

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# variables
CPU=30
EXTEND=200

# usage
USAGETXT=\
"
$0 <macs2 peak bed file> <genome fasta> <meme db> <output directory>
"

# test 
#isEnvVarSet $UPSCb

# load modules
module load bioinfo-tools BEDTools samtools MEME

# Arguments: bed-file genome-fasta database-path outdir
# The input should be: 1 - bed-file, 2 - genome-fasta, 3 - database-path, 4 - outdir

if [ ! -f $(which bedtools) ]; then
    abort "bedtools is not available"
fi

if [ ! -f $(which meme-chip) ]; then
    abort "meme-chip is not available"
fi

if [ ! -f $1 ]; then
    abort "The first argument needs to be a bed file that contains a peak."
fi

if [ ! -f $2 ]; then
    abort "The second argument needs to be the Populus genome in .fasta format."
fi

if [ ! -f $2.fai ]; then
   abort "The genome should have been indexed in .fai format. Use samtools faidx to do so."
fi

if [ ! -f $3 ]; then
    abort "The third argument needs to be the path to the database file containing known binding motifs."
fi

if [ ! -d $4 ]; then
    mkdir $4
fi

# designate temporary files
SIZE=$(tempfile)
FLANK=$(tempfile)
GENOME=$2
FASTA=$(tempfile)


# create the GENOME file. Takes the full Potra GENOME .fai file and extracts the SIZE of the chromosomes.
# This enables defining the FLANKs (eg. doesn't make up imaginary sequences near the ends of the chromosomes)
cut -f1,2 ${GENOME}.fai > $SIZE

# for a bed file. This adds -200 [1bp PEAK] +200 FLANKs to the peak and creates a .bed file of the result.
bedtools flank -l $EXTEND -r $EXTEND -i $1 -g $SIZE > $FLANK
# or flank -b 200

bedtools getfasta -fi $GENOME -bed $FLANK -fo $FASTA

# Finally, this runs meme-chip and makes output of the motifs
meme-chip -db $3 $FASTA -oc $4

