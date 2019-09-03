#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e
set -x

## modules
module load bioinfo-tools pasa BEDTools

## a usage function
usage(){
    echo >&2 \
"
    Usage: $0 [options] <database name> <pasa dir>

    Options:
             -C provide a comparison ID to extract, default to latest
             -A exports all PASA assemblies including lncRNAs
             -I add to include a an intergenic filter gff to be created
             -g add a gff reference

    Note: this prints to the stdout! 

" 
    exit 1
}

## GLOBAL OPTIONS
COMPID=
DOALL=0
IG=0
gff="$2/genome.gff3"

while getopts C:AIg: option
 do
     case "$option" in
   C) COMPID="-C $OPTARG";;
   A) DOALL=1;;
    I) IG=1;;
    g) gff="$OPTARG";;
         \?) ## unknown flag
             usage;;
   esac
done
shift `expr $OPTIND - 1`


## we get one name and a dir as input
if [ $# != 2 ]; then
    echo "This function takes a database name and an existing pasa dir as argument"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory"
    usage
    if [[ ! -f genome || ! -f transcripts.fasta.clean ]]; then
    echo "This does not look like a pasa results directory"
    usage
    fi
fi

## cd in the out dir
cd $2

## submit
if [ $DOALL -eq 0 ]; then
    $PASAHOME/scripts/dump_valid_annot_updates.dbi -M $1 -V -R -g ../genome $COMPID
else 
    ## to retrieve all sequences
    $PASAHOME/scripts/dump_valid_annot_updates_modified_v1.dbi -M $1 -g ../genome $COMPID -A > pasa-dump-full.gff3
    if [ $IG -eq 0 ]; then
        awk '{if ( $3 == "gene" ){ print $0 }}' pasa-dump-full.gff3 > full_output_gene_locus.gff3
        bedtools subtract -A -a full_output_gene_locus.gff3 -b $gff > full_non_overlapping_gene_locus.gff3
        bedtools intersect -wa -f 1 -u -a pasa-dump-full.gff3 -b full_non_overlapping_gene_locus.gff3 > pasa-dump-intergenic.gff3

        ## clean up
        rm full_non_overlapping_gene_locus.gff3

        ##  TODO: Pasa validator and a script to remove sub features without parent should be added

    fi
fi 
