#! /bin/bash -l
#SBATCH -p core -n 1
#SBATCH -t 0-04:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e

module load bioinfo-tools
module load blast/2.2.27+
module load samtools/0.1.19

if [ $# != 2 ]; then
    echo "The argument should be the database name to download and the output directory"
    exit 1
fi

db=
case "$1" in
    nt) db="nt";;
    nr) db="nr";;
esac

if [ -z $db ]; then
    echo "The argument should be one of nt or nr"
    exit 1
fi  

if [ ! -d $2 ]; then
    echo "The  directory name you provided does not exist"
    exit 1
fi

##
echo "Setting up"

## copy the exec to proper dir
cp $SLURM_SUBMIT_DIR/../../../src/perl/update_blastdb.pl $2

## cd to that dir
cd $2

## download
echo "Downloading"
$2/update_blastdb.pl --passive --force $1

## unpacking
echo "Unpacking"
find . -name "*.tar.gz" -exec tar -zxf "{}" \;
## TODO to go parallel
## find . -name "*.tar.gz" -print0 | xargs -P 16 -I {} -0  tar -zxf {}

## cleaning
echo "Cleaning"
rm *.tar.gz

## extracting the fasta seq
echo "Extracting the fasta"
blastdbcmd -db $1 -entry all -out $1.fa
samtools faidx $1.fa
awk 'BEGIN{FS=" "};{print $2,$1}' $1.fa.fai > $1.lengths.txt

## add a flag
echo "Flagging"
touch ${1}-Updated-`date +%Y%m%d`.flag

echo "Done"

