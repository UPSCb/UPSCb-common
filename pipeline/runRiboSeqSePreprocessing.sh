#!/bin/bash
set -eu

### ========================================================
# Preprocessing script for RNA-Seq data.
# THIS SCRIPT IS NOT TO BE RUN THROUGH SBATCH, USE BASH!
### ========================================================
VERSION="0.0.4"

### ========================================================
## pretty print
### ========================================================
underline=`tput smul`
nounderline=`tput rmul`
bold=`tput bold`
normal=`tput sgr0`

### ========================================================
## functions
### ========================================================

## usage
usage() {
    echo "usage: bash `basename $0` [OPTIONS] <proj> <mail> <fastq> <outdir>

Run the RiboSeq preprocessing pipeline, i.e. FastQC, Trimmomatic, Sortmerna, Bowtie and kallisto. 
You throw in a fastq file and it spits out a pair of BAM files. Sweet!

${bold}ARGUMENTS:${normal}
    proj    project that this should be run as
    mail    e-mail address where notifications will be sent
    fastq   single-end fastq file
    outdir  directory to save everything in

${bold}OPTIONS:${normal}
    -h        show this help message and exit
    -d        dry-run. Check for the software dependencies, output the command lines and exit
    -D        Debug. On Uppmax use a devel node to try out changes during development or program updates, it's highly recommended to not run more than one step in this mode
    -s n      step at which to start (see ${underline}STEPS${nounderline})
    -e n      step at which to end (see ${underline}STEPS${nounderline})
    -p dir    the provisional directory (default to /tmp or $SNIC_TMP) otherwise
    -b dir    path to Bowtie index (required if Bowtie is included in pipeline)
    -m        the memory requirement for performing the alignment, in GB; i.e. 128 (default) allocates 128GB
    -r        the reference databases for sortmerna
    -c        clipping arguments paased to trimmomatic. Overrides defaults in runTrimmomatic.sh. Default in this script: ILLUMINACLIP:\"$TRIMMOMATIC_HOME/adapters/NEB-universal-adapter.fa\":1:15:10
    -T        trimming arguments passed to trimmomatic. Overrides defaults in runTrimmomatic.sh. Default in this script (we do not trim based on quality) : MINLEN:16
    -f fasta  the transcript fasta sequence file for kallisto
    -k inx    the Kallisto index file (use the index that has a kmer size of 15)
    -M        fragment length mean for kallisto
    -S        fragment length sd for kallisto
    -t        library is s${underline}t${nounderline}randed (illumina protocol, second strand cDNA using dUTP)
    -a        library is s${underline}t${nounderline}randed (non illumina protocol e.g. first strand cDNA using dUTP)

    
${bold}STEPS:${normal}
    The steps of this script are as follows:

    1) FastQC
    2) SortMeRNA
    3) FastQC
    4) Trimmomatic
    5) FastQC
    6) Bowtie
    7) Kallisto
    
    Tentative extensions
    8) MultiQC
    9) Ribowaltz

${bold}NOTES:${normal}
    * This script should not be run through sbatch, just use bash.
    * If the output directory already exists, content may be overwritten.
    * All tools needed must be in your PATH. You will be notified if they aren't.
    * If the SORTMERNA and UPSCb environmental variables don't exist, the
      script will guess them. Set them yourself just to be safe.
    * If sbatch is not available, then the script does not submit the jobs
      but rather print out the commands.
" 1>&2
    exit 1
}

## check if a tool is  present and is executable
toolCheck() {
    tool=`which $1 2>/dev/null`
    if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
	echo 0
    else
	echo 1
    fi
}

## cleanup on failed test
cleanup() {
    #jobs=$@
    if [ $# -gt 0 ]; then
        echo "Canceling already started jobs: $@" 1>&2
        scancel $@
    fi
    exit 3
}

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
    time=""

    while getopts "J:e:o:d:Dm:t:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            D) debug="-p devel -t 1:00:00";;
	          m) memory=$OPTARG ;;
	          t) time=$OPTARG ;;
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
    ${memory:+"--mem $memory"} \
    ${dependency:+"-d" "$dependency"} \
    ${time:+"-t $time"} \
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
    time=""

    while getopts "J:e:o:d:Dm:t:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            D) debug="";;
            m) memory=$OPTARG ;;
            t) time=$OPTARG ;;
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

### ========================================================
## main
### ========================================================

## set the default vars
dryrun=0
debug=0
pstart=1
pend=7
mem=
memArg="--mem"
trimmomatic_clipping=ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/NEB-universal-adapter.fa:1:15:10
trimmomatic_options="MINLEN:16"
bowtie_index=
ilm_stranded=0
non_ilm_stranded=0
kallisto_index=
kallisto_fasta=
fragment_length_mean=
fragment_length_sd=
sortmerna_db=

# Parse the options
OPTIND=1
while getopts "ab:c:Dde:f:hk:M:m:p:r:S:s:tT:" opt; do
    case "$opt" in
        h) usage;;
        a) non_ilm_stranded=1;;
        b) bowtie_index=$OPTARG ;;
        c) trimmomatic_clipping=$OPTARG ;;
	    d) dryrun=1;;
        D) debug=1;;
        s) pstart=$OPTARG ;;
        e) pend=$OPTARG ;;
        f) kallisto_fasta=$OPTARG ;;
        k) kallisto_index=$OPTARG ;;
	m) mem="${OPTARG}GB";;
	p) SNIC_TMP=$OPTARG ;;
	    r) sortmerna_db="-d $OPTARG";;
	    t) ilm_stranded=1 ;;
	    T) trimmomatic_options=$OPTARG ;;
	    M) fragment_length_mean=$OPTARG ;;
	    S) fragment_length_sd=$OPTARG ;;
        ?) usage;;
    esac
done
shift $((OPTIND - 1))

# check for exclusive options
if [ $ilm_stranded -eq 1 ] && [ $non_ilm_stranded -eq 1 ]; then
    echo "ERROR: the -a and -t option are mutually exclusive. Choose either or."
    usage
fi

# check the number of arguments
if [ "$#" !=  "4" ]; then
    echo "ERROR: this script expects 4 arguments"
    usage
fi

## Do basic checks
[[ $pstart =~ ^[0-9]+$ ]] || {
    echo "ERROR: '$pstart' is not a valid start value" 1>&2
    usage
}
[[ $pend =~ ^[0-9]+$ ]] || {
    echo "ERROR: '$pend' is not a valid end value" 1>&2
    usage
}
[[ $pend -lt $pstart  ]] && {
    echo "ERROR: you can't finish before you start" 1>&2
    usage
}


echo "### ========================================
# UPSC RiboSeq pre-processing pipeline v$VERSION
### ========================================
# Checking the environment for the necessary tools"

## trimmotatic is java based and we have the jar 
## in the git repos, so we only need to check
## for java
## the toolList list all the necessary tools
## the toolArray (starting at 1) link the tool(s) to its respective step(s)
## starArray and htseqArray are there to simulate a nested array
toolList=(fastqc sortmerna java bowtie2 samtools kallisto)
bowtieArray=([0]=3 [1]=4)
kallistoArray=([0]=4 [1]=5)
toolArray=([1]=0 [2]=1 [3]=0 [4]=2 [5]=0 [6]=${bowtieArray[*]} [7]=${kallistoArray[*]})


## a global var to stop if we miss tools
## but only after having checked them all
status=0
for ((i=$pstart;i<=$pend;i++)); do
    ## the nested arrays have a length of (2 x number of element) -1
	tot=$(((${#toolArray[$i]} + 1) / 2))
	for ((j=0;j<$tot;j++)); do
	    echo "# Checking step $i ${toolList[${toolArray[$i]:$j:$((j + 1))}]}"
	    if [ `toolCheck ${toolList[${toolArray[$i]:$j:$((j + 1))}]}` -eq 1 ]; then
		echo >&2 "ERROR: Please install the missing tool: ${toolList[${toolArray[$i]:$j:$((j + 1))}]}."
		((status+=1))
	    fi
	done
done
if [ $status -gt 0 ]; then
    exit 1;
fi

echo "### ========================================
# Preparation started on `date`
# Going from step $pstart to $pend
# Provided parameters are:
# bowtie index: $bowtie_index
# Kallisto index: $kallisto_index"

# Check if bowtie is included and if the reference is set
if [ $pstart -le 6 ] && [ $pend -ge 6 ]; then
    if [ -z $bowtie_index ]; then
        echo >&2 "ERROR: You are running bowtie but have not given an index (-b)"
        usage
    elif [ ! -d $(dirname bowtie_index) ]; then
        echo >&2 "ERROR: Could not find the bowtie index directory"
        usage
    fi
fi

# Check if kallisto is included and if the files are provided
if [ $pstart -le 7 ] && [ $pend -ge 7 ]; then
    if [ -z $kallisto_index ]; then
        echo >&2 "ERROR: You are running kallisto but have not given an index (-k)"
        usage
    elif [ ! -f $kallisto_index ]; then
        echo >&2 "ERROR: Could not find the kallisto index"
        usage
    fi
    if [ -z $kallisto_fasta ]; then
        echo >&2 "ERROR: You are running kallisto but have not given the fasta reference (-f)"
        usage
    elif [ ! -f $kallisto_fasta ]; then
        echo >&2 "ERROR: Could not find the kallisto fasta reference"
        usage
    fi
fi

# This variable holds the absolute path of this script, i.e. the
# repository's pipeline directory.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## abort on any error
set -eu

## get some vars
proj=$1
mail=$2

## check some more
if [ ! -f $3 ]; then
    echo "ERROR: fastq is not a file: '$3'" 1>&2
    usage
    usage
fi

## Just to make sure we avoid path issues
# fastq1=$(readlink -f "$3")
fastq1="$3"

# Sample name to use for output
s_prefix=${fastq1%.f*q.gz}
sname=`basename $s_prefix`

if [ ! -d `dirname $4` ]; then
    echo "ERROR: could not find parent directory for output directory" 1>&2
    usage
    usage
fi

## Set up the directory structure
## in case outdir is .
outdir=`readlink -f $4`
[[ ! -d $outdir ]] && mkdir $outdir

## Use the directory name as a job identifier
JOBNAME=`basename $outdir`

## FastQC
fastqc_raw="$outdir/fastqc/raw"
[[ ! -d $fastqc_raw ]] && mkdir -p $fastqc_raw
fastqc_sortmerna="$outdir/fastqc/sortmerna"
[[ ! -d $fastqc_sortmerna ]] && mkdir $fastqc_sortmerna
fastqc_trimmomatic="$outdir/fastqc/trimmomatic"
[[ ! -d $fastqc_trimmomatic ]] && mkdir $fastqc_trimmomatic

## SortmeRNA
sortmerna="$outdir/sortmerna"
[[ ! -d $sortmerna ]] && mkdir $sortmerna

## Trimmomatic
trimmomatic="$outdir/trimmomatic"
[[ ! -d $trimmomatic ]] && mkdir $trimmomatic

## Bowtie
bowtie="$outdir/bowtie2"
[[ ! -d $bowtie ]] && mkdir $bowtie

## Kallisto
kallisto="$outdir/kallisto"
[[ ! -d $kallisto ]] && mkdir $kallisto

## Export some variables
: "${UPSCb:=$PIPELINE_DIR/..}"
#[[ -z $UPSCb ]] && export UPSCb=$PIPELINE_DIR/..

## final setup
## check for sbatch, return 1 if no sbatch, 0 otherwise
CMD=sbatch
if [ `toolCheck $CMD` -eq 1 ]; then
    CMD=bash
fi

## setup the tmp dir
#[[ -z $SNIC_TMP ]] && export SNIC_TMP=/tmp
# this makes the above set -u friendly
# From https://stackoverflow.com/questions/11362250/in-bash-how-do-i-test-if-a-variable-is-defined-in-u-mode
: "${SNIC_TMP:=/tmp}"

echo "### ========================================"

## Check if we dry run
if [ $dryrun -eq 1 ]; then
    echo "# This is a dry-run! Printing command lines only"
fi
echo "# Using $CMD as a job submitter"

### --------------------------------------------------------
## prepare the job array
debug_var=
if [ $debug -eq 1 ]; then
        echo "# Running in debug mode, devel node will be used on uppmax, it's highly recommended to not run more than one step in this mode"
        debug_var="-D"
fi

## Job ID array
JOBIDS=()

## Job CMD array
JOBCMDS=()

## Internal field separator
## We need to change the input field separator 
## to be line feed and not space for some of the
## array, but this might fail command execution
SPACEIFS=$IFS
LFIFS=$'\n'
IFS=$LFIFS

## TODO WE NEED A single prepare and run function
## that contains the logic
## Run FastQC
if [ $pstart -le 1 ]; then
    echo "# Preparing step 1"

    JOBCMDS+=(`prepare_$CMD \
        -e $fastqc_raw/${sname}_1_fastqc.err \
        -o $fastqc_raw/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.raw1 \
        runFastQC.sh $fastqc_raw $fastq1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

fi

## Run Sortmerna. Depends on a successful run of FastQC
if [ $pstart -le 2 ] && [ $pend -ge 2 ]; then
    echo "# Preparing step 2"

    # because we can have no unset variable, we first try a parameter expansion
    # it set the variable to empty if it does not exist
    SORTMERNADIR=${SORTMERNADIR:-}
    if [ -z $SORTMERNADIR ]; then
        echo "ERROR: could not find the sortmerna data in $SORTMERNADIR" 1>&2
        usage
    fi
    
    dep=
    if [ $pstart -lt 2 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    fi

    JOBCMDS+=(`prepare_$CMD \
        $debug_var \
        -e $sortmerna/${sname}_sortmerna.err \
        -o $sortmerna/${sname}_sortmerna.out \
        $dep \
        -J ${sname}.RNAseq.SortMeRNA \
	-t 2-00:00:00 \
        runSortmerna.sh $sortmerna_db -P 11,6,1 -s 1 -a -u $sortmerna $SNIC_TMP $fastq1`)
    if [ "$CMD" == "bash" ]; then
	JOBCMDS+=("export SORTMERNADIR=$SORTMERNADIR;$sortmerna_id")
    else
	JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
    fi
fi

fastq_sort_1="$sortmerna/${sname}_sortmerna.fq.gz"

# Run FastQC. Depends on a successful run of SortMeRNA
if [ $pstart -le 3 ] && [ $pend -ge 3 ]; then
    echo "# Preparing step 3"

    dep=
    if [ $pstart -lt 3 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    elif [ ! -f $fastq_sort_1 ] && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-file could not be found"
        cleanup
    fi

    JOBCMDS+=(`prepare_$CMD \
        -e $fastqc_sortmerna/${sname}_1_fastqc.err \
        -o $fastqc_sortmerna/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.SortMeRNA1 \
        $dep runFastQC.sh $fastqc_sortmerna $fastq_sort_1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

## Run trimmomatic. Depends on a successful run of SortMeRNA
if [ $pstart -le 4 ] && [ $pend -ge 4 ]; then
    echo "# Preparing step 4"
    
    dep=
    if [ $pstart -le 2  ] && [ "$CMD" != "bash" ]; then
        # SortMeRNA has to finish, if it was started
	    # hence we check for a pstart lower or equal to 3
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-2]}"
    elif [ ! -f $fastq_sort_1 ] && [ "$CMD" != "bash" ]; then
        # If we start here, make sure the files exist
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    JOBCMDS+=(`prepare_$CMD \
        $debug_var \
        -e $trimmomatic/${sname}_trimmomatic.err \
        -o $trimmomatic/${sname}_trimmomatic.log \
        -J ${sname}.RNAseq.Trimmomatic \
        $dep \
        runTrimmomatic.sh -t -s -c $trimmomatic_clipping $fastq_sort_1 $trimmomatic $trimmomatic_options`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

fi

## Trimmed fastq file paths
fastq_trimmed_1="$trimmomatic/${sname}_sortmerna_trimmomatic.fq.gz"

# Run FastQC. Depends on a successful run of Trimmomatic
if [ $pstart -le 5 ] && [ $pend -ge 5 ]; then
    echo "# Preparing step 5"

    dep=
    if [ $pstart -lt 4 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    elif [ ! -f $fastq_trimmed_1 ] && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: Trimmed FASTQ-files could not be found"
        cleanup
    fi
    JOBCMDS+=(`prepare_$CMD \
        -e $fastqc_trimmomatic/${sname}_1_fastqc.err \
        -o $fastqc_trimmomatic/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.Trimmomatic1 \
        $dep runFastQC.sh $fastqc_trimmomatic $fastq_trimmed_1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

fi

# Run Bowtie. Depends on a successful run of Trimmomatic.
if [ $pstart -le 6 ] && [ $pend -ge 6 ]; then
    echo "# Preparing step 6"

    dep=
    ## Trimmomatic needs finishing if it was started
    ## hence, we check if step 5 was started  
    if [ $pstart -le 4 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-2]}"
    elif [ ! -f $fastq_trimmed_1 ] && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: Trimmed FASTQ-files could not be found"
        cleanup
    fi

    JOBCMDS+=(`prepare_$CMD \
            $debug_var \
            -e $bowtie/${sname}_bowtie.err \
            -o $bowtie/${sname}_bowtie.out \
            -J ${sname}.RNAseq.bowtie \
            ${mem:+"-m" "$mem"} \
            -t 2-00:00:00 \
            $dep runBowtie2.sh -s $bowtie_index $fastq_trimmed_1 $bowtie $SNIC_TMP`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

# Run Kallisto. Depends on a successful run of trimmomatic
if [ $pstart -le 7 ] && [ $pend -ge 7 ]; then
    echo "# Preparing step 7"

    dep=
    if [ $pstart -le 4 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-3]}"
    elif ([ ! -f $fastq_trimmed_1 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: Trimmed FASTQ-file could not be found"
        cleanup
    fi
    
    kallisto_strand_option=
    if [ $ilm_stranded -eq 0 ] && [ $non_ilm_stranded -eq 0 ]; then
	    kallisto_strand_option="-u" 
    elif [ $ilm_stranded -eq 1 ]; then
      kallisto_strand_option="-r"   
    elif [ $non_ilm_stranded -eq 1 ]; then
	    kallisto_strand_option="-F"
    fi
  
    single_end_option="-s -M $fragment_length_mean -S $fragment_length_sd"
    JOBCMDS+=(`prepare_$CMD \
        $debug_var \
        -e $kallisto/${sname}_kallisto.err \
        -o $kallisto/${sname}_kallisto.out \
        -J ${sname}.RNAseq.kallisto \
        $dep runKallisto.sh -p $single_end_option $kallisto_strand_option $fastq_trimmed_1 $kallisto_index $kallisto_fasta $kallisto`)
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
	echo "${JOBCMDS[*]}";;
esac
