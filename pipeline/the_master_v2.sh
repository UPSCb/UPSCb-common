#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -A b2010064

##
##  TODO  continue from automatic .success
##

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
script_dir="/home/davidsu/UPSCb/pipeline"
validate=false
## Input data
source_dir=$1"/"$5
out_dir=$2"/"$5
fq_out=$2
sname=$3
fread=$3"_1.fq.gz"
rread=$3"_2.fq.gz"
genomeFolder=$4
gff=$genomeFolder"/Ptrichocarpa_210_gene_exons.gff3"
sample=$5
proj=$6
trimming=$7
force=false
countFile=$sample"_count_lostreads.txt"

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


##  Run fastQC on raw files
if [ ! -d "$fq_out/fastQC" ]; then
    mkdir $fq_out/fastQC
    mkdir $fq_out/fastQC/raw
    mkdir $fq_out/fastQC/trimmomatic
    mkdir $fq_out/fastQC/sortmerna
fi

###  fastq validation
if $validate; then
    echo "Fastq validation" >> $logfile
    sh $script_dir/runFastQValidator.sh $source_dir/$fread
    sh $script_dir/runFastQValidator.sh $source_dir/$rread
    echo "Validation OK" >> $logfile

    echo "Run fastqc" >> $logfile

    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/raw $source_dir/$fread 
    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/raw $source_dir/$rread
###                    ###

fi


###   Sortmerna  ####
echo "Run sortmerna"  >> $logfile
if [ ! -d "$out_dir/tmp" ]; then
    mkdir $out_dir/tmp
fi
if [ ! -d $out_dir/sortmerna ]; then
    mkdir $out_dir/sortmerna
fi

tmp_dir=$out_dir/tmp

sh $script_dir/runSortmerna.sh $out_dir/sortmerna $tmp_dir $source_dir/$fread $source_dir/$rread

#../python/repair2.0.py $out $in1 $in2

fread=${fread%_1.fq.gz}_sortmerna_1.fq.gz
rread=${rread%_2.fq.gz}_sortmerna_2.fq.gz

source_dir=$out_dir/sortmerna

## fastq validation
if $validate; then
    echo "Fastq validation" >> $logfile
    sh $script_dir/runFastQValidator.sh $source_dir/$fread
    sh $script_dir/runFastQValidator.sh $source_dir/$rread
    echo "Validation OK" >> $logfile
    ##  Run fastQC on files
    echo "Run fastqc" >> $logfile

    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/sortmerna $source_dir/$fread 
    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/sortmerna $source_dir/$rread

fi


######################

### Trimmomatic ####

echo "Run trimmomatic" >> $logfile

if [ ! -d "$out_dir/trimmomatic" ]; then
    mkdir $out_dir/trimmomatic
fi

# Usage: $0 <fwd fastq file> <rev fastq file> <output dir> [trimming option]"


sh $script_dir/runTrimmomatic.sh $source_dir/$fread $source_dir/$rread $out_dir/trimmomatic $trimming
 

fread=${fread%_1.fq.gz}_trimmomatic_1.fq.gz
rread=${rread%_2.fq.gz}_trimmomatic_2.fq.gz

source_dir=$out_dir/trimmomatic

## fastq validation
if $validate; then
    echo "Fastq validation" >> $logfile
    sh $script_dir/runFastQValidator.sh $source_dir/$fread
    sh $script_dir/runFastQValidator.sh $source_dir/$rread
    echo "Validation OK" >> $logfile

    ##  Run fastQC on files
    echo "run FastQC" >> $logfile
    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/trimmomatic $source_dir/$fread 
    sbatch -A $proj $script_dir/runFastQC.sh $fq_out/fastQC/trimmomatic $source_dir/$rread

fi


####################
 
#########
##
##   Start alignment.
##
#########


####    STAR    ####

## usage  $prog -a $fread -b $rread -o $dirOut -m 11000 -g $gff -t $genemodel
echo "run STAR" >> $logfile

if [ ! -d "$out_dir/STAR" ]; then
    echo "here"
    mkdir $out_dir/STAR
fi

minintron=11000

sh $script_dir/runSTAR.sh $source_dir/$fread $source_dir/$rread $genomeFolder/STAR $gff -m $minintron -o "$out_dir/STAR"

alignment=${fread%_1.fq.gz}_STAR.bam

echo "Alignment complete" >> $logfile
####################

##### Cufflinks #####
echo "run Cufflinks" >> $logfile
sh $script_dir/runCufflinks.sh $out_dir $genomeFolder/genome.fa $alignment
echo "cufflinks complete" >> $logfile
####################

echo "Done!" >> $logfile