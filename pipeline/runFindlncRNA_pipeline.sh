#!/bin/bash
set -e

### ========================================================
# Postprocessing script for RNA-Seq data, find long non coding RNA sequences and annotate.
# THIS SCRIPT IS NOT TO BE RUN THROUGH SBATCH, USE BASH!
### ========================================================
VERSION="0.0.1"


### ========================================================
## pretty print
### ========================================================
underline=`tput smul`
nounderline=`tput rmul`
bold=`tput bold`
normal=`tput sgr0`

usage(){
	echo '''

runFindlncRNA_pipeline.sh [OPTIONS] <proj> <mail> <pasafolder> <pasadatabase> <genome> <gff3>

${bold}OBSERVE!!!!!${normal} ----  Under development!

${bold}STEPS:${normal}
    The steps of this script are as follows:

    1) Extract PASA EST sequences
    2) Filter intergenic genes (optional, false by default) (will be run with 1)
    3) get FASTA sequences, both prot and nucleotide.
    4) frameDP
	5) filter sequences
		a) full length sequences
		b) 40% CDS
	6) Step has to be done in R, filter by expression

${bold}Arguments:${normal}
	required:
		<proj>			project Name
		<mail>			user mail
		<dir>			Path to the PASA folder 
		<pasaDB>		Name of the PASA database
		<genome>		Path to the species genome fasta
	optional:
		-C 				Specify a database version for PASA default (latest)
		-i 				Add to look for only intergenic novel genes.
		-d 				dry run
		-c 				specify frameDP CFG file, default cfg in P.tricocarpa ref
		-t 				Train dataset for frameDP, in case no previous training has 
						been done on the species. takes path to the species fasta file T also needs to be set
        -T              give a path to the frameDP traning set reference folder (Required to be set together with -t)
		-s 				Step to start with
		-e 				Step to end with
'''
}

#usage

[[ -z $UPSCb ]] && export UPSCb=$PIPELINE_DIR/..

# This variable holds the absolute path of this script, i.e. the
# repository's pipeline directory.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

### --------------------------------------------------------
## sbatch related 
run_sbatch_usage() {
    echo "usage: run_sbatch OPTIONS <batch script> [<script param> ...]

    Run a batch script and echo the job id

    OPTIONS
      -e  stderr (required)
      -o  stdout (required)
      -J  job name
      -d  dependency
      -D  debug mode
      -m memory" 1>&2
    exit 1
}

prepare_sbatch() {
    # Start a batch script and echo the job id
    if [ $# -lt 3 ]; then
        echo "ERROR: wrong number of arguments"
        run_sbatch_usage
    fi

    OPTIND=1

    log_path=""
    out_path=""
    dependency=""
    memory=""
    debug=""

    while getopts "J:e:o:d:Dm:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            D) debug="-p devel -t 1:00:00";;
            m) memory=$OPTARG ;;
            ?) run_sbatch_usage;;
        esac
    done

    shift $((OPTIND-1))

    script=$1
    shift

    if [ -z $jobname ]; then
        jobname="${sname}.RNAPreproc.${script}"
    fi

    # Check that the output file paths are given
    if [ -z $log_path ] || [ -z $out_path ]; then
        echo "ERROR: output file paths are not given"
        run_sbatch_usage
    usage
    fi
    # Check that the output file directories exist
    if [ ! -d `dirname $log_path` ] || [ ! -d `dirname $out_path` ]; then
        echo "ERROR: stderr and stdout paths must exist" 1>&2
        run_sbatch_usage
    usage
    fi

    # sbatch_options=
    # if [ ! -z $memory ]; then
    # sbatch_options=" -C $memory"
    # fi

    # if [ ! -z $dependency ]; then
    # sbatch_options="${sbatch_options} -d $dependency"
    # fi

    # echo the command
    echo "sbatch -A" $proj \
    "-J" $jobname "--mail-user" $mail \
    "-e" $log_path "-o" $out_path \
    $debug \
    ${memory:+"-C" "$memory"} \
    ${dependency:+"-d" "$dependency"} \
    $PIPELINE_DIR/$script $@
    
}

## TODO warn about the random!!! 
run_sbatch() {
    if [ $dryrun -eq 1 ]; then
    echo $RANDOM
    else
    IFS=$SPACEIFS
    sbatch_echo=`$1`
    if [ $? -ne 0 ]; then
            echo "ERROR: submission failed" 1>&2
        echo "$1" 1>&2
            cleanup ${JOBIDS[*]}
    fi
    IFS=$LFIFS
    echo ${sbatch_echo//[^0-9]/}
    fi
}

### --------------------------------------------------------
## bash commands
prepare_bash() {

    # echo the job
    # fetch but ignore the arguments
    # 
    OPTIND=1
    log_path=""
    out_path=""
    dependency=""
    memory=""
    debug=""

    while getopts "J:e:o:d:D:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            D) debug="";;
            ?) run_sbatch_usage;;
        esac
    done
    shift $((OPTIND-1))

    # the script
    script=$1
    shift

    echo "bash $PIPELINE_DIR/$script $@;"
}

run_bash(){
    echo $RANDOM
}


####  First step create run command for the PASA dump script, remember the -A option.

CFG="/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/cfg/FrameDP_new.cfg"

## Get optional arguments
OPTIND=1
INTERGENIC=""
dbversion=""
dryrun=0
pstart=1
pend=4
species_fasta=""
frameDP_ref=""

while getopts "C:c:t:s:e:id" opt; do
	case "$opt" in
		C) dbversion=" -C $OPTARG";;
		c) CFG="$OPTARG";;
		t) species_fasta="$OPTARG";;
        T) frameDP_ref="$OPTARG";;
		s) pstart="$OPTARG";;
		e) pend="$OPTARG";;
		i) INTERGENIC=" -I";;
		d) dryrun=1;;
	esac
done
shift $((OPTIND-1))

## Get required arguments
proj=$1
mail=$2

folder=$3
pasaDB=$4
genome=$5
gff=$6

## Other defaults
sname="pasa_novel_genes"
fname=$folder"/novel_genes"
if [ ! -d $fname ]; then
    mkdir $fname
fi


toolCheck() {
    tool=`which $1 2>/dev/null`
    if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
	echo 0
    else
	echo 1
    fi
}

CMD=sbatch
if [ `toolCheck $CMD` -eq 1 ]; then
    CMD=bash
fi
## command array
JOBCMD=();
JOBIDS=();

####  Second step, have to be run optionally 

if [ $pstart -le 1 ]; then
    echo "# Preparing step 1 (and 2) (Extract EST sequences)"
    opt="$INTERGENIC"
    opt+="$dbversion"
    JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_dump_EST_gff.err \
        -o $fname/${sname}_dump_EST_gff.out \
        -J ${sname}.Extract_EST \
        runPasaValidAnnotationUpdate.sh -A$opt -g $gff $pasaDB $fname`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

if [ $pstart -le 2 ] && [ $pend -ge 2 ]; then
    echo "# Preparing step 3 (extract FASTA seq)"
    if [ "$INTERGENIC" != "" ]; then
    	pasa_gff="pasa-dump-intergenic"
    else
    	pasa_gff="pasa-dump-full"
    fi
    JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_extract_fasta.err \
        -o $fname/${sname}_extract_fasta.out \
        -J ${sname}.Extract \
        runGetPasaFasta.sh $fname/${pasa_gff}.gff3 $genome $fname/${pasa_gff}.fasta`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

if [ $pstart -le 3 ] && [ $pend -ge 3 ]; then
    echo "# Preparing step 4 (frameDP) and NCBI blast"
    if [ $species_fasta != "" ]; then
    	echo "Create a training set if it doesnt exist"
    	JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_frameDP_train.err \
        -o $fname/${sname}_frameDP_train.out \
        -J ${sname}.frameDP \
        runFrameDP.sh $species_fasta $folder/frameDP $CFG $frameDP_ref --only_train`)
    	JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
    fi
    echo "running data"
    JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_frameDP.err \
        -o $fname/${sname}_frameDP.out \
        -J ${sname}.frameDP \
        runFrameDP.sh $folder/${pasa_gff}.fasta $folder/frameDP $CFG --notrain`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
    echo "running BLAST ncbi"
        JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_BLAST_ncbi.err \
        -o $fname/${sname}_frameDP_train.out \
        -J ${sname}.frameDP \
        runBlastPlusSplit.sh $species_fasta $folder/frameDP $CFG `)
        JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

if [ $pstart -le 4 ] && [ $pend -ge 4 ]; then
	##  Filter steps including merging all result files from frameDP, 
	##	finding fl sequences, filter CDS shorter than 40% of the mRNA etc.
    echo "# Preparing step 4 (Filtering steps)"
    JOBCMDS+=(`prepare_$CMD \
        -e $fname/${sname}_1_validate.err \
        -o $fname/${sname}_1_validate.out \
        -J ${sname}.PostFiltering \
        runPythonPostFiltering.sh `)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

echo "### ========================================
# Preparation done on `date`"

## done
case $CMD in
    sbatch) 
	if [ $dryrun -eq 1 ];then
	    echo "${JOBCMDS[*]}"
	else
	    echo "Successfully submitted ${#JOBIDS[@]} jobs for ${sname}: ${JOBIDS[@]}"
	fi;;
    bash)
	## because IFS was set to '\n', we need to use ${JOBCMDS[*]} instead of ${JOBCMDS[@]}
	echo "${JOBCMDS[@]}";;
esac