#!/bin/bash -l
#SBATCH -p node 
#SBATCH -n 8
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=ALL

set -e #stop on error

if [ "$1" == "help" ]; then
	echo "Script:  Master script with all steps to complete the alignment pipeline."
	echo "sortmerna -> trimmomatic -> STAR -> cufflinks -> cuffmerge+cuffcompare"
	echo "Usage: <base_dir_in> <base_dir_out> <readname> <genome_folder> <sample(tree)> <proj>"
	echo "ex. the_master.sh /proj/b2010064/sequence_data/raw/K1 /proj/b2010064/nobackup/STARpipeline/K1 K1-01 Ptrichocarpa_210 K1"
	echo "In the base_dir_in reads will be stored after sortmerna+trimmomatic, and STAR.bam and parsedCuffdiff and cuffmerge data."
	echo "Other data will be temporary stored in nobackup"
fi

##  General paths
data_dir="/proj/$6/nobackup/data"
script_dir="~/UPSCb/pipeline"

## Input data
#source_dir=$1
out_dir=$2
sname=$3
fread=$3"_1.fq.gz"
rread=$3"_2.fq.gz"
genomeFolder=$4
sample=$5
trimming=$6
source_dir=$1"/"$5

if [ "$trimming" == "" ]; then
    trimming="SLIDINGWINDOW:5:30 MINLEN:50"
fi

if [ ! -d "$out_dir" ]; then
    mkdir $out_dir
fi

## create master logfile
if [ ! -d "$out_dir/logs" ]; then
    mkdir $out_dir/logs
fi
touch $out_dir"/logs/"$fread".log"
logfile=$out_dir"/logs/"$fread".log"

echo "Star running sample: $sample " > $logfile

echo ""
echo "Running options" >> $logfile
echo "    source_dir: "$source_dir >> $logfile
echo "    out_dir: "$out_dir >> $logfile
echo "    Reads: "$fread" "$rread >> $logfile
echo "    genomeFolder: "$genomeFolder >> $logfile
echo "    trimming opt: "$trimming >> $logfile
echo ""


###  fastq validation

echo "Fastq validation" >> $logfile
echo sh $script_dir/runFastQValidator.sh $source_dir/$fread
echo sh $script_dir/runFastQValidator.sh $source_dir/$rread
echo "Validation OK" >> $logfile

##  Run fastQC on raw files
if [ ! -d "$out_dir/fastQC" ]; then
    mkdir $out_dir/fastQC
fi

echo "Run fastqc" >> $logfile

echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/raw $source_dir/$fread
echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/raw $source_dir/$rread
###                    ###

###   Sortmerna  ####
echo "Run sortmerna"  >> $logfile
if [ ! -d $out_dir/tmp ]; then
    mkdir $out_dir/tmp
fi
if [ ! -d $out_dir/sortmerna ]; then
    mkdir $out_dir/sortmerna
fi

tmp_dir=$out_dir/tmp

echo $script_dir/runSortmerna.sh $out_dir/sortmerna $tmp_dir $fread $rread

fread=${fread%_1.fq.gz}_sortmerna_1.fq.gz
rread=${rread%_2.fq.gz}_sortmerna_2.fq.gz

source_dir=$out_dir/sortmerna

## fastq validation
echo "Fastq validation" >> $logfile
echo sh $script_dir/runFastQValidator.sh $source_dir/$fread
echo sh $script_dir/runFastQValidator.sh $source_dir/$rread
echo "Validation OK" >> $logfile
##  Run fastQC on files

echo "Run fastqc" >> $logfile

echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/sortmerna $source_dir/$fread
echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/sortmerna $source_dir/$rread

######################

### Trimmomatic ####

echo "Run trimmomatic" >> $logfile

if [ ! -d $out_dir/trimmomatic ]; then
    mkdir $out_dir/trimmomatic
fi

# Usage: $0 <fwd fastq file> <rev fastq file> <output dir> [trimming option]"

echo $script_dir/runTrimmomatic.sh $source_dir/$fread $source_dir/$rread $out_dir/trimmomatic $trimming

fread=${fread%_1.fq.gz}_trimmomatic_1.fq.gz
rread=${rread%_2.fq.gz}_trimmomatic_2.fq.gz

source_dir=$out_dir/trimmomatic

## fastq validation
echo "Fastq validation" >> $logfile
echo sh $script_dir/runFastQValidator.sh $source_dir/$fread
echo sh $script_dir/runFastQValidator.sh $source_dir/$rread
echo "Validation OK" >> $logfile
##  Run fastQC on files
echo "run FastQC" >> $logfile
echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/trimmomatic $source_dir/$fread
echo sbatch $script_dir/runFastQC.sh $out_dir/fastQC/trimmomatic $source_dir/$rread

####################

#########
##
##   Start alignment.
##
#########


####    STAR    ####

## usage  $prog -a $fread -b $rread -o $dirOut -m 11000 -g $gff -t $genemodel
echo "run STAR" >> $logfile

if [ ! -d $out_dir/STAR ]; then
    mkdir $out_dir/STAR
fi

minintron=11000

echo sh $script_dir/runSTAR.sh -a $source_dir/$fread -b $source_dir/$rread $out_dir/STAR -m $minintron -t $genomeFolder/STAR

alignment=${fread%_1.fq.gz}_STAR.bam

echo "Alignment complete" >> $logfile
####################

##### Cufflinks #####
echo "run Cufflinks" >> $logfile
echo sh $script_dir/runCufflinks.sh $out_dir $genomeFolder/genome.fa $alignment
echo "cufflinks complete" >> $logfile
####################

echo "Done!" >> $logfile