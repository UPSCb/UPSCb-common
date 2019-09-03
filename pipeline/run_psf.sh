#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -ex

# load the modules
module load bioinfo-tools psf

# usage function
usage(){
echo >&2 \
"
	Usage: psf <list(.cfg)> <seq_prot(.fa)> <-psf_cfg:psf.cfg> <-o:pmap.cfg> <-O:prot_mapM.CFG> <-search_here:seq_nuc.fa> <-no_echo> 

		name of the executable file ./psf + ...

		list(.cfg)	             - path to the list file
		seq_prot.fa              - path to the multiFASTA-file with protein sequences, without gaps. Headers can include additional information in Softberry 						   AbInitio or FGENESH++ format. Here IPI or NR database could be given on input.
		-psf_cfg:psf.cfg	       - path to the psf_cfg file (psf configuration file)
		-o:pmap.cfg              - path to configuration file with parameters of the general alignment algorithm
		-O:prot_mapM.CFG         - path to configuration file with the options of general protein-on-DNA mapping algorithm
		-search_here:seq_nuc.fa  - path nucleotide FASTA-file with a single genomic sequence (without gaps).
		-no_echo
"
	exit 1
}

# check the arguments

if [ $# != 7 ]; then
    echo "This function requires 7 arguments."
    usage
fi


# run the command
psf $1 $2 $3 $4 $5 $6 $7 


