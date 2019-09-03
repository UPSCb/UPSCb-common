#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 62
#SBATCH --mem 500GB
#SBATCH	--mail-user=ALL
#SBATCH -w watson

set -e
#set -x

echo "Usage: runLoadPasaGenome.sh PASAconfig.config genome.fasta genome.gff"

module load bioinfo-tools pasa

config=$1
fasta=$2
gff=$3
trinity=$4
cufflinks=$5
s=25

## Validate file

#$PASAHOME/misc_utilities/pasa_gff3_validator.pl $gff

#if [ $? == 0 ]; then
#	echo $gff" GFF validation ok"
#else
#	echo "GFF validation failed please verify that your gff file is correct"
#	exit
#fi

$PASAHOME/scripts/Launch_PASA_pipeline.pl -c $config -s 25 -R -g $fasta -t $trinity --ALIGNERS blat,gmap --CPU 31 --cufflinks_gtf $cufflinks

## create a new database
#./scripts/create_mysql_cdnaassembly_db.dbi -c config -S $PASAHOME/schema/cdna_alignment_mysqlschema

# load current genome
#$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi -c $config -g $fasta -P $gff

echo "Done!"
#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 62
#SBATCH --mem 500GB
#SBATCH	--mail-user=ALL
#SBATCH -w watson

PASAHOME="/mnt/picea/Modules/apps/bioinfo/pasa/r20140417"

echo "Usage: runLoadPasaGenome.sh PASAconfig.config genome.fasta genome.gff"

module load bioinfo-tools
module load pasa/r20140417 blat/36 trinity/r20140717

config=$1
fasta=$2
gff=$3
trinity=$4
cufflinks=$5

#echo $1 $2 $3 $4 $5

#exit

## Validate file


#$PASAHOME/misc_utilities/pasa_gff3_validator.pl $gff

#if [ $? == 0 ]; then
#	echo $gff" GFF validation ok"
#else
#	echo "GFF validation failed please verify that your gff file is correct"
#	exit
#fi

$PASAHOME/scripts/Launch_PASA_pipeline.pl -c $config -s 2 -R -g $fasta -t $trinity --ALIGNERS blat,gmap --CPU 31 --cufflinks_gtf $cufflinks

## create a new database
#./scripts/create_mysql_cdnaassembly_db.dbi -c config -S $PASAHOME/schema/cdna_alignment_mysqlschema

# load current genome
#$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi -c $config -g $fasta -P $gff

echo "Done!"
