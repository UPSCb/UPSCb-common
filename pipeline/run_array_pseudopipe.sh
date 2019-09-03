#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x

# load the modules
module load bioinfo-tools 
module load Pseudopipe

# usage function
usage(){
echo >&2 \
"
	Usage: $0 <output dir> <masked dna file path> <input dna file path> <input pep file path> <exon mask file path> <PBS?>

		- <output dir>: the directory where the output files will be created

		- <input dna dir>: contains a file named dna_rm.fa, which is entire repeat masked dna from that species, 
		and a list file for all unmasked dna divided into different chromosomes in FASTA format;

		- <input pep dir>: contains a FASTA file for all the proteins in the species;

		- <exon mask dir>: contains a list of files named \"chr_exLocs\", \"chr2_exLocs\", etc. to specify exons coordinates, one for each chromosome. 
		Only thing matters for these files are their third and fourth columns, which should be start and end coordinates of exons.

		- <PBS> : 0 if, 1 if.

"
	exit 1
}

# check values

if [ $# != 6 ]; then
    echo "This function requires 6 arguments."
    usage
fi

if [ ! -d $1 ]; then
	echo "The output directory : $1 does not exist"
	usage
fi

if [ ! -f $2 ]; then
	echo "the masked dna sequence file $2 does not exist"
	usage
fi

if [ "${2##*.}" != ".fa" ]; then
  echo "The file needs to be a fasta file, i.e. have a .fa extension"
  usage
fi

if [ "${3##*.}" != ".fa" ]; then
  echo "The file needs to be a fasta file, i.e. have a .fa extension"
  usage
fi

if [ ! -f $4 ]; then
	echo "the protein sequence file $4 does not exist"
	usage
fi

if [ "${4##*.}" != ".fa" ]; then
  echo "The file needs to be a fasta file, i.e. have a .fa extension"
  usage
fi

if [ "${5##*.}" != "_exLocs" ]; then
  echo "The file needs to end by \"_exLocs\""
  usage
fi

if [ ! $6 ]; then
	echo "The argument: $6 does not exist"
	usage
fi


# run the command
#pseudopipe.sh $1 $2 $3.$SLURM_ARRAY_TASK_ID $4 $5 $6
pseudopipe.sh $1 $2 $3 $4 $5 $6

