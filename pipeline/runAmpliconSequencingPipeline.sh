#!/bin/bash
set -ex

### ========================================================
# Preprocessing script for RNA-Seq data.
# THIS SCRIPT IS NOT TO BE RUN THROUGH SBATCH, USE BASH!
### ========================================================
# DEFAULTS
VERSION="0.0.9"

ACCOUNT=
ADAPTER=
ADAPTER16S=$UPSCb/data/adapter16S.fa
ADAPTERITS=$UPSCb/data/adapterITS.fa
ADAPTERBOTH=$UPSCb/data/adapterboth.fa
CLOSED=0
EMAIL=
SRNA="16S"
CPU=16
MEM=2
TRIM="SLIDINGWINDOW:20:20 MINLEN:100"
PARAMETER=
REF=

### ========================================================
## functions
### ========================================================
# usage function
usage(){
  echo >&2 \
  "
  Usage: bash $(basename $0) [options] <R1.fastq.gz> <R2.fastq.gz> <I1.fastq.gz> <I2.fastq.gz> <map file> <out dir> 
  
  Options: 
    -a adapter file (if not 16S / ITS default)
    -A SLURM account
    -c closed OTU picking (default to open)
    -h print this message and exits
    -e SLURM email
    -m SLURM --mem requirement (the given value will be multiplied by the number of CPUs (-t) and is in GB)
    -p parameter file for OTU picking
    -r reference file for OTU picking
    -s '16S' or 'ITS' (default to 16S)
    -t SLURM -n : maximum number of cpu per node to use per job (job array are used for further parallelisation)
    -T trimmomatic trimming option - default to $TRIM
  
  Notes:
    The UPSCb environment variable needs to be defined and point to the root of your UPSCb checkout from Git.
    
  Details: The script runs
    1. FastQC (raw)
    2. DeML
    3. FastQC (DeML)
    4. Trimmomatic
    5. FastQC (Trimmomatic)
    6. Flash
    7. FastQC (Flash)
    8. MultiQC
    9. ...
    
    
  "
  exit 1
}

abort(){
  echo >&2 $1
  usage
}

# sbatch
run_sbatch_usage() {
    echo "usage: run_sbatch OPTIONS <batch script> [<script param> ...]

    Run a batch script and echo the job id

    OPTIONS
      -a  job array specification
      -e  stderr (required)
      -o  stdout (required)
      -J  job name
      -n  number of CPU
      -m  memory to use
      -d  dependency" 1>&2
}

cleanup() {
    #jobs=$@
    if [ $# -gt 0 ]; then
        echo "Canceling already started jobs: $@" 1>&2
        scancel $@
    fi
    exit 3
}

run_sbatch() {
    # Start a batch script and echo the job id
    if [ $# -lt 3 ]; then
        run_sbatch_usage
        exit 1
    fi

    OPTIND=1

    log_path=""
    out_path=""
    dependency=""
    array=""
    node="2"
    mem=$(expr $node "*" 2)
    

    while getopts "J:a:e:m:n:o:d:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            a) array=$OPTARG;;
            e) log_path=$OPTARG ;;
            m) mem=$OPTARG ;;
            n) node=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            ?) run_sbatch_usage; exit 1 ;;
        esac
    done

    shift $((OPTIND-1))

    script=$1
    shift

    if [ -z $jobname ]; then
        jobname="${sname}.$(basename $script)"
    fi
    
    # Check that the output file paths are given
    if [ -z $log_path ] || [ -z $out_path ]; then
        run_sbatch_usage
        exit 1
    fi
    # Check that the output file directories exist
    if [ ! -d `dirname $log_path` ] || [ ! -d `dirname $out_path` ]; then
        echo "ERROR: stderr and stdout paths must exist" 1>&2
        run_sbatch_usage
        exit 1
    fi

    sbatch_options=
    if [ ! -z $dependency ]; then
        sbatch_options="-d $dependency"
    fi

    if [ ! -z $array ]; then
        sbatch_options="$sbatch_options -a $array"
    fi

    sbatch_echo=`sbatch -A "$ACCOUNT" \
                -J "$jobname" \
                --mail-user "$EMAIL" \
                --mail-type=ALL \
                --mem=${mem}GB \
                -n $node \
                -e "$log_path" \
                -o "$out_path" \
                $sbatch_options \
                $script $@`

    if [ $? -ne 0 ]; then
        echo "ERROR: submission failed" 1>&2
        cleanup ${JOBIDS[*]}
    fi

    echo ${sbatch_echo//[^0-9]/}
}

### ========================================================
## main
### ========================================================

# check that the UPSCb env. var. exists
if [ -z $UPSCb ] || [ ! -d $UPSCb ]; then
  abort "The UPSCb env. var. needs to be defined"
  usage
fi

# load the modules
module load bioinfo-tools FastQC deML Trimmomatic flash/1.2.11 fastx multiqc # Qiime blast/2.2.26 python/2.7.8

# process the options
while getopts "a:A:ce:hm:p:r:s:t:T:" opt; do
  case $opt in
	a) ADAPTER=$OPTARG;;
	A) ACCOUNT=$OPTARG;;
	c) CLOSED=1;;
	e) EMAIL=$OPTARG;;
	h) usage;;
	m) MEM=$OPTARG;;
	p) PARAMETER=$OPTARG;;
	r) REF=$OPTARG;;
	s) SRNA=$OPTARG;;
  t) CPU=$OPTARG;;
  T) TRIM=$OPTARG;;
  *) abort "Unknown option: $opt";;
  esac
done
shift `expr $OPTIND - 1`

# check the options
# -----------------
case $SRNA in
  16S) 
    if [ -z $ADAPTER ]; then
        ADAPTER=$ADAPTER16S
    fi
  ;;
  ITS) 
    if [ -z $ADAPTER ]; then
        ADAPTER=$ADAPTERITS
    fi
  ;;
  Both)
      if [ -z $ADAPTER ]; then
	  ADAPTER=$ADAPTERBOTH
      fi
  ;;
  *) abort "Only valid values for -s are '16S' and 'ITS'";;
esac

if [ ! -f $ADAPTER ]; then
  abort "The adapter file $ADAPTER does not exist"
fi

if [ ! -z $PARAMETER ] && [ ! -f $PARAMETER ]; then
  abort "The parameter file $PARAMETER does not exist"
fi

if [ ! -z $REF ] && [ ! -f $REF ]; then
  abort "The reference file $REF does not exist"
fi

# check the arguments
# -------------------
if [ $# -ne 6 ]; then
  abort "The script expect 6 arguments"
fi
     
if [ ! -f $1 ]; then
 abort "The forward fastq file: $1 does not exist"
fi
F=$1
     
if [ ! -f $2 ]; then
 abort "The reverse fastq file: $2 does not exist"
fi
R=$2
    
if [ ! -f $3 ]; then
  abort "The forward index file: $3 does not exist"
fi
I1=$3

if [ ! -f $4 ]; then
 abort "The reverse index file: $4 does not exist"
fi
I2=$4

if [ ! -f $5 ]; then
 abort "The mapping file: $5 does not exist"
fi
M=$5

if [ ! -d $6 ]; then
  abort "The output directory: $6 does not exist"
fi
OUT=$(readlink -f $6/)

# declare vars
# ------------
declare -A toolsArray
declare -A cpuArray
declare -A memArray
toolsArray=(
  ["asArray"]="$UPSCb/pipeline/runAsArray.sh"
  ["FastQC"]="$UPSCb/pipeline/runFastQC.sh"
  ["DeML"]="$UPSCb/pipeline/runDeML.sh"
  ["Trimmomatic"]="$UPSCb/pipeline/runTrimmomatic.sh"
  ["Flash"]="$UPSCb/pipeline/runFLASH.sh"
  ["MultiQC"]="$UPSCb/pipeline/runMultiQC.sh"
)

# our cluster defaults to 2 CPU
cpuArray=(
  ["FastQC"]="2"
  ["DeML"]="2"
  ["Trimmomatic"]="$CPU"
  ["Flash"]="$CPU"
  ["MultiQC"]="2"
)

memArray=(
  ["FastQC"]="$(expr ${cpuArray["FastQC"]} "*" 2)"
  ["DeML"]="$(expr ${cpuArray["DeML"]} "*" 2)"
  ["Trimmomatic"]="$(expr ${cpuArray["Trimmomatic"]} "*" $MEM)"
  ["Flash"]="$(expr ${cpuArray["Flash"]} "*" $MEM)"
  ["MultiQC"]="$(expr ${cpuArray["MultiQC"]} "*" 2)"
)

# Setup the out dir
# -------------------
mkdir -p $OUT/Biom
mkdir -p $OUT/DeML
mkdir -p $OUT/FastQC/raw
mkdir -p $OUT/FastQC/DeML
mkdir -p $OUT/FastQC/Flash
mkdir -p $OUT/FastQC/Trimmomatic
mkdir -p $OUT/Flash
mkdir -p $OUT/MultiQC
mkdir -p $OUT/Trimmomatic

sname=$(basename ${F//.fastq.gz/})
smpNumber=$(expr $(cat $M | wc -l) - 1)
jobNumber=$(expr $smpNumber - 1)

# define the job array step lists
step03_list=$OUT/FastQC/DeML/list.txt
step04_list=$OUT/Trimmomatic/list.txt
step05_list=$OUT/FastQC/Trimmomatic/list.txt
step06_list=$OUT/Flash/list.txt
step07_list=$OUT/FastQC/Flash/list.txt

# remove previous run
rm -f $step03_list
rm -f $step04_list
rm -f $step05_list
rm -f $step06_list
rm -f $step07_list

# fill them in
for n in $(awk '{if(NR >1){print $3}}' $M); do

  # step 3
  echo "$OUT/FastQC/DeML $OUT/DeML/demultiplex_${n}_r1.fq.gz $OUT/DeML/demultiplex_${n}_r1.fq.gz" >> $step03_list  

  # step 4
  ln -sf $OUT/DeML/demultiplex_${n}_r1.fq.gz $OUT/Trimmomatic/demultiplex_${n}_1.fq.gz
  ln -sf $OUT/DeML/demultiplex_${n}_r2.fq.gz $OUT/Trimmomatic/demultiplex_${n}_2.fq.gz
  echo "-c "ILLUMINACLIP:${ADAPTER}:2:30:10" -p ${cpuArray["Trimmomatic"]} $OUT/Trimmomatic/demultiplex_${n}_1.fq.gz $OUT/Trimmomatic/demultiplex_${n}_2.fq.gz" >> $step04_list

  # step 5
  echo "$OUT/FastQC/Trimmomatic $OUT/Trimmomatic/demultiplex_${n}_trimmomatic_1.fq.gz $OUT/Trimmomatic/demultiplex_${n}_trimmomatic_2.fq.gz" >> $step05_list
  
  # step 6
  echo "-t ${cpuArray["Flash"]} -c -o $OUT/Trimmomatic/demultiplex_${n}_trimmomatic_1.fq.gz $OUT/Trimmomatic/demultiplex_${n}_trimmomatic_2.fq.gz" >> $step06_list

  # step 7
  echo "$OUT/FastQC/Flash $OUT/Flash/demultiplex_${n}_trimmomatic_flash.fq.gz" >> $step07_list

done

# Queue the jobs
# --------------
JOBIDS=()

# Run raw FastQC
step01=$OUT/FastQC/raw/${sname}_fastqc.zip
if [ ! -f $step01 ]; then
    step01_id=$(run_sbatch \
        -e $OUT/FastQC/raw/${sname}.err \
        -o $OUT/FastQC/raw/${sname}.out \
        -n ${cpuArray["FastQC"]} -m ${memArray["FastQC"]} \
        ${toolsArray["FastQC"]} $OUT/FastQC/raw $F $R)
    JOBIDS+=($step01_id)
fi

# Run DeML
step02=$(find $OUT/DeML -name "demultiplex*_i1.fq.gz" | wc -l)
if [ $step02 -ne $smpNumber ]; then
    step02_id=$(run_sbatch \
        -e $OUT/DeML/${sname}.err \
        -o $OUT/DeML/${sname}.out \
        -n ${cpuArray["DeML"]} -m ${memArray["DeML"]} \
        ${toolsArray["DeML"]} $F $R $I1 $I2 $M $OUT/DeML/demultiplex)
    JOBIDS+=($step02_id)
fi

# Run DeML FastQC (job array)
step03=$(find $OUT/FastQC/DeML -name "demultiplex*_r1_fastqc.zip" | wc -l)
if [ $step03 -ne $smpNumber ]; then
    dep=
    if [ ! -z $step02_id ]; then
      dep="-d aftercorr:$step02_id"
    fi
    step03_id=$(run_sbatch \
        -e $OUT/FastQC/DeML/${sname}%a.err \
        -o $OUT/FastQC/DeML/${sname}%a.out \
        $dep \
        -a "0-$jobNumber" \
        -n ${cpuArray["FastQC"]} -m ${memArray["FastQC"]} \
        ${toolsArray["asArray"]} ${toolsArray["FastQC"]} $step03_list)
    JOBIDS+=($step03_id)
fi

# Run Trimmomatic
step04=$(find $OUT/Trimmomatic -name "demultiplex*_trimmomatic_1.fq.gz" | wc -l)
if [ $step04 -ne $smpNumber ]; then
    dep=
    if [ ! -z $step02_id ]; then
      dep="-d afterok:$step02_id"
    fi
    step04_id=$(run_sbatch \
        -e $OUT/Trimmomatic/${sname}%a.err \
        -o $OUT/Trimmomatic/${sname}%a.out \
        $dep \
        -a "0-$jobNumber" \
        -n ${cpuArray["Trimmomatic"]} -m ${memArray["Trimmomatic"]} \
        ${toolsArray["asArray"]} ${toolsArray["Trimmomatic"]} $step04_list $OUT/Trimmomatic $TRIM)
    JOBIDS+=($step04_id)
fi

# Run Trimmomatic FastQC
step05=$(find $OUT/FastQC/Trimmomatic -name "demultiplex*_trimmomatic_1_fastqc.zip" | wc -l)
if [ $step05 -ne $smpNumber ]; then
    dep=
    if [ ! -z $step04_id ]; then
      dep="-d aftercorr:$step04_id"
    fi
    step05_id=$(run_sbatch \
        -e $OUT/FastQC/Trimmomatic/${sname}%a.err \
        -o $OUT/FastQC/Trimmomatic/${sname}%a.out \
        $dep \
        -a "0-$jobNumber" \
        -n ${cpuArray["FastQC"]} -m ${memArray["FastQC"]} \
        ${toolsArray["asArray"]} ${toolsArray["FastQC"]} $step05_list)
    JOBIDS+=($step05_id)
fi

# Run Flash
step06=$(find $OUT/Flash -name "demultiplex*_trimmomatic_flash.fa" | wc -l)
if [ $step06 -ne $smpNumber ]; then
    dep=
    if [ ! -z $step04_id ]; then
      dep="-d aftercorr:$step04_id"
    fi
    step06_id=$(run_sbatch \
        -e $OUT/Flash/${sname}%a.err \
        -o $OUT/Flash/${sname}%a.out \
        $dep \
        -a "0-$jobNumber" \
        -n ${cpuArray["Flash"]} -m ${memArray["Flash"]} \
        ${toolsArray["asArray"]} ${toolsArray["Flash"]} $step06_list $OUT/Flash -M 251 -O)
    JOBIDS+=($step06_id)
fi

# Run FastQC Flash
step07=$(find $OUT/FastQC/Flash -name "demultiplex*_trimmomatic_flash_fastqc.zip" | wc -l)
if [ $step07 -ne $smpNumber ]; then
    dep=
    if [ ! -z $step06_id ]; then
      dep="-d aftercorr:$step06_id"
    fi
    step07_id=$(run_sbatch \
        -e $OUT/FastQC/Flash/${sname}%a.err \
        -o $OUT/FastQC/Flash/${sname}%a.out \
        $dep \
        -a "0-$jobNumber" \
        -n ${cpuArray["FastQC"]} -m ${memArray["FastQC"]} \
        ${toolsArray["asArray"]} ${toolsArray["FastQC"]} $step07_list)
    JOBIDS+=($step07_id)
fi

# Run multiqc
if [ ! -f $OUT/MultiQC/multiqc_report.html ]; then
    dep=
    if [ ! -z $step07_id ]; then
      dep="-d aftercorr:$step07_id"
    fi
    step08_id=$(run_sbatch \
        -e $OUT/MultiQC/${sname}.err \
        -o $OUT/MultiQC/${sname}.out \
        $dep \
        -n ${cpuArray["MultiQC"]} -m ${memArray["MultiQC"]} \
        ${toolsArray["MultiQC"]} $OUT $OUT/MultiQC)
    JOBIDS+=($step08_id)
fi

     # 
     # #ITS picking with a closed reference
     # if ["$8"=="closed"]; then
     # mkdir -p $7/OTUclosed
     # for f in $(find $7/fasta -type f -name \
     #            "*.fq.gz.fa"); 
     # do
     # smpldir=$7/OTU/$(basename ${f//.fq.gz.fa/})
     # mkdir $smpldir
     # pick_closed_reference_otus.py -i $f -f \
     # -r $9 \
     # -o $smpldir --assign_taxonomy -p $10 parameterfile.txt; #--parallel -O 12;
     # done
     # fi
     # 
     # if ["$8"=="open"]; then
     # mkdir -p $7/OTUopen
     # for f in $(find $7/fasta -type f -name \
     #            "*.fq.gz.fa"); 
     # do
     # smpldir=$7/OTU/$(basename ${f//.fq.gz.fa/})
     # mkdir $smpldir
     # pick_closed_reference_otus.py -i $f -f \
     # -r $9 \
     # -o $smpldir --assign_taxonomy -p $10 parameterfile.txt; #--parallel -O 12;
     # done
     # fi
     # 
     # #Converts otu table from to a biom file
     # for f in $(find $7/OTUs -type f -name "otu_table.biom");
     # do 
     # pattern=$(basename $(dirname $f));
     # mkdir -p $7/biom/$pattern;
     # biom convert -i $f -o $7/biom/$pattern/table.from_biom_w_taxonomy.txt --to-tsv --table-type "OTU table" --header-key taxonomy;
     # done
     # 
