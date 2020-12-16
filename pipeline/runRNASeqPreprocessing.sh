#!/bin/bash
set -eu

### ========================================================
# Preprocessing script for RNA-Seq data.
# THIS SCRIPT IS NOT TO BE RUN THROUGH SBATCH, USE BASH!
### ========================================================
VERSION="0.3.0"

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
    echo "usage: bash `basename $0` [OPTIONS] <proj> <mail> <fastq1> <fastq2> <outdir>

Run the RNA-Seq preprocessing pipeline, i.e. FastQC, Trimmomatic, Sortmerna
and Salmon. You throw in a pair of fastq files and it spits out an sf file. Sweet!

${bold}ARGUMENTS:${normal}
    proj    project that this should be run as
    mail    e-mail address where notifications will be sent
    fastq1  forward fastq file
    fastq2  reverse fastq file
    outdir  directory to save everything in

${bold}OPTIONS:${normal}
    -h        show this help message and exit
    -d        dry-run. Check for the software dependencies, output the command lines and exit
    -D        Debug. On Uppmax use a devel node to try out changes during development or program updates, it's highly recommended to not run more than one step in this mode
    -s n      step at which to start (see ${underline}STEPS${nounderline})
    -e n      step at which to end (see ${underline}STEPS${nounderline})
    -E        do not skip the step 1,8-10 (skipped by default)
    -f fasta  the transcript fasta sequence file for kallisto
    -F fasta  the genome fasta file (required if STAR is included in the pipeline)
    -g dir    path to STAR reference to use (required if STAR is included
              in pipeline)
    -k        the STAR genome is available in shared memory (the STAR option ${underline}--genomeLoad LoadAndKeep${nounderline})
    -l        the BAM sorting memory limit for STAR in bytes (default: 10000000000)
    -K inx    the Kallisto index file
    -G gtf    gene model GTF file for STAR
    -H gff    gff3 file for ${underline}H${nounderline}TSeq
    -i        IDATTR in GFF3 file to report counts (default: 'Parent')
    -m        the memory requirement for performing the alignment (in GB; e.g. enter 128 for 128GB)
    -p        fastq data is phred64 encoded
    -t        library is s${underline}t${nounderline}randed (currently only relevant for HTSeq, for the illumina protocol, second strand cDNA using dUTP)
    -a        library is s${underline}t${nounderline}randed (currently only relevant for HTSeq, for the non illumina protocol e.g. first strand cDNA using dUTP)
    -S        Salmon index
    -T        trimming arguments passed to trimmomatic. Overrides defaults in runTrimmomatic.sh
    -I        the max intron length for STAR
    -x        the sortmerna index directory

${bold}STEPS:${normal}
    The steps of this script are as follows:

    1) (fastQValidator*)
    2) FastQC
    3) SortMeRNA
    4) FastQC
    5) Trimmomatic
    6) FastQC
    7) Salmon
    8) Kallisto*
    9) STAR*
    10) (HTSeq-count*)
    
    The steps marked with * are optional; only accessible if -E is set.
    The steps in () are tools not supported in kogia

${bold}NOTES:${normal}
    * This script should not be run through sbatch, just use bash.
    * If the output directory already exists, content may be overwritten.
    * All tools needed must be in your PATH. You will be notified if they aren't.
    * If the SORTMERNA and UPSCb environmental variables don't exist, the
      script will guess them. Set them yourself just to be safe.
    * If sbatch is not available, then the script does not submit the jobs
      but rather print out the commands.
    * The -a and -t option are mutually exclusive
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
    ${memory:+"$memArg $memory"} \
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

    while getopts "J:e:o:d:Dm:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            D) debug="";;
        m) memory=$OPTARG ;;
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
# This variable holds the absolute path of this script, i.e. the
# repository's pipeline directory.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## set the default vars
dryrun=0
debug=0
pstart=1
pend=9
star_ref=
star_gtf=
star_runner_options=
star_options=
htseq_gff=
ilm_stranded=0
non_ilm_stranded=0
idattr="Parent"
mem=
memArg="--mem"
bam_memory=10000000000
phred_value=
trimmomatic_options=
kallisto_index=
kallisto_fasta=
star_intron_max=
genome_fasta=
salmon_index=
sortmerna_inx=
extended=0
kogia=${PIPELINE_DIR}/../kogia/scripts

# Parse the options
OPTIND=1
while getopts "adDe:Ef:F:g:G:hH:i:I:kK:l:m:ps:S:tT:x:" opt; do
    case "$opt" in
        a) non_ilm_stranded=1;;
	      d) dryrun=1;;
        D) debug=1;;
        e) pend=$OPTARG ;;
        E) extended=1;;
        f) kallisto_fasta=$OPTARG ;;
        F) genome_fasta=$OPTARG ;;
        g) star_ref=$OPTARG ;;
        G) star_gtf=$OPTARG ;;
        h) usage;;
        H) htseq_gff=$OPTARG ;;
        i) idattr=$OPTARG ;;
        I) star_intron_max=$OPTARG ;;
        k) star_options="-- --genomeLoad LoadAndKeep";;
        K) kallisto_index=$OPTARG ;;
        l) bam_memory=$OPTARG ;;
        m) mem="${OPTARG}GB";;
        p) phred_value="-q";;
        s) pstart=$OPTARG ;;
        S) salmon_index=$OPTARG ;;
        t) ilm_stranded=1 ;;
		    T) trimmomatic_options=$OPTARG ;;
		    x) sortmerna_inx="-i $OPTARG";;
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
if [ "$#" !=  "5" ]; then
    echo "ERROR: this script expects 5 arguments"
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

[[ $extended -eq 0 ] && [[ $pstart -eq 1 ] || [ $pstart -gt 7 ]]] && {
  echo "ERROR: you cannot run step 1,8-10 unless you set -E" 1>&2
  usage
}


echo "### ========================================
# UPSC pre-RNA-Seq pipeline v$VERSION
### ========================================
# Checking the environment for the necessary tools"

## trimmotatic is java based and we have the jar 
## in the git repos, so we only need to check
## for java
## the toolList list all the necessary tools
## the toolArray (starting at 1) link the tool(s) to its respective step(s)
## starArray and htseqArray are there to simulate a nested array
toolList=(fastQValidator fastqc sortmerna trimmomatic salmon kallisto samtools star htseq-count)
htseqArray=([0]=6 [1]=8)
kallistoArray=([0]=5 [1]=6)
toolArray=([1]=0 [2]=1 [3]=2 [4]=1 [5]=3 [6]=1 [7]=4 [8]=${kallistoArray[*]} [9]=7 [10]=${htseqArray[*]})

## a global var to stop if we miss tools
## but only after having checked them all
status=0
for ((i=$pstart;i<=$pend;i++)); do
    ## the nested arrays have a length of (2 x number of element) -1
	tot=$(((${#toolArray[$i]} + 1) / 2))
	for ((j=0;j<$tot;j++)); do
	    echo "# Checking step $i ${toolList[${toolArray[$i]:$j:$((j + 1))}]}"
	    if [ $(toolCheck $kogia/${toolList[${toolArray[$i]:$j:$((j + 1))}]}) -eq 1 ] && [ $(toolCheck ${toolList[${toolArray[$i]:$j:$((j + 1))}]}) -eq 1 ]; then
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
# Salmone index: $salmon_index
# Kallisto index: $kallisto_index
# Kallisto fasta: $kallisto_fasta
# STAR genome: $star_ref 
# STAR gtf: $star_gtf
# HTSeq gff3: $htseq_gff"

# Check if salmon is included and then if the reference is set
if [ $pstart -le 7 ] && [ $pend -ge 7 ]; then
   if [ -z $salmon_index ]; then
        echo >&2 "ERROR: You are running salmon but have not given an index, set the -S option"
        usage
    elif [ ! -d $salmon_index ]; then
        echo >&2 "ERROR: Could not find the salmon index directory"
        usage
    fi
fi

if [ $extended -eq 1 ]; then

  # Check if kallisto is included and then if the reference is set
  if [ $pstart -le 8 ] && [ $pend -ge 8 ]; then
     if [ -z $kallisto_index ]; then
          echo >&2 "ERROR: You are running kallisto but have not given an index, set the -K option"
          usage
      elif [ ! -f $kallisto_index ]; then
          echo >&2 "ERROR: Could not find the kallisto index file"
          usage
      fi
      if [ -z $kallisto_fasta ]; then
          echo >&2 "ERROR: You are running kallisto but have not given a fasta reference, set the -f option"
          usage
      elif [ ! -f $kallisto_fasta ]; then
          echo >&2 "ERROR: Could not find the kallisto fasta file"
          usage
      fi 
  fi

  # Check if STAR is included and then if the reference is set
  if [ $pstart -le 9 ] && [ $pend -ge 9 ]; then
    if [ -z $star_ref ]; then
        echo >&2 "ERROR: You are running STAR but have not given a STAR reference, set the -g option"
        usage
    elif [ ! -f $star_ref/Genome ]; then
        echo >&2 "ERROR: Could not find STAR reference"
        usage
    fi
    if [ ! -z $star_gtf ] && [ ! -f $star_gtf ]; then
        echo >&2 "ERROR: Could not find gtf file: '$star_gtf'"
        usage
    fi
    
    if [ ! -z $genome_fasta ] && [ ! -f $genome_fasta ]; then
        echo >&2 "ERROR: Could not find the genome fasta file: '$genome_fasta', set the -F option"
        usage
    fi
    
    # Check that STAR will get a GTF file
    if [ ! -z $star_gtf ] && [ -f $star_gtf ]; then
        if [ ${star_gtf##*.} != "gtf" ]; then
          echo >&2 "ERROR: You should use a gene model GTF file for STAR"
          usage
        fi
    fi
  fi

  # Check that HTSeq will get a GFF3 file (if it's included in the pipeline)
  if [ $pstart -le 10 ] && [ $pend -ge 10 ]; then
    if [ -z $htseq_gff ] && [ -z $star_gtf ]; then
        echo >&2 "ERROR: HTSeq needs a GFF3 file"
        usage
    elif [ -z $htseq_gff ]; then
        htseq_gff=$star_gtf
    elif [ ! -z $htseq_gff ] && [ ! -f $htseq_gff ]; then
        echo >&2 "ERROR: Could not find gff file: '$htseq_gff'"
        usage
    fi
  fi
fi

## get some vars
proj=$1
mail=$2

## check some more
if [ ! -f $3 ]; then
    echo "ERROR: fastq1 is not a file: '$3'" 1>&2
    usage
fi

if [ ! -f $4 ]; then
    echo "ERROR: fastq2 is not a file: '$4'" 1>&2
    usage
fi

## Just to make sure we avoid path issues
# fastq1=$(readlink -f "$3")
# fastq2=$(readlink -f "$4")
fastq1="$3"
fastq2="$4"

# Sample name to use for output
s_prefix=${fastq1%_1.f*q.gz}
sname=$(basename $s_prefix)

if [ ! -d $(dirname $5) ]; then
    echo "ERROR: could not find parent directory for output directory" 1>&2
    usage
fi

## Set up the directory structure
## in case outdir is .
outdir=$(readlink -f $5)
[[ ! -d $outdir ]] && mkdir $outdir

## Use the directory name as a job identifier
JOBNAME=$(basename $outdir)

dirList=(fastqc/raw fastqc/raw sortmerna fastqc/sortmerna trimmomatic fastqc/trimmomatic salmon kallisto star htseq-count)

for ((i=$pstart;i<=$pend;i++)); do
  echo "# Creating directory for step $i: $outdir/$dirList[$((i -1))]" 
	mkdir -p $outdir/$dirList[$((i -1))]
done

## final setup
## check for sbatch, return 1 if no sbatch, 0 otherwise
CMD=sbatch
if [ `toolCheck $CMD` -eq 1 ]; then
    CMD=bash
fi

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

## Run fastQValidator
step=0
if [ $pstart -le (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_1_validate.err \
        -o $dirList[$step]/${sname}_1_validate.out \
        -J ${sname}.RNAseq.FastQValidate1 \
        runFastQValidator.sh $fastq1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_2_validate.err \
        -o $dirList[$step]/${sname}_2_validate.out \
        -J ${sname}.RNAseq.FastQCValidate2 \
        runFastQValidator.sh $fastq2`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

## Run FastQC. Depends on a successful run of fastQValidator
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep1=
    dep2=
    if [ $pstart -lt (($step+1)) ] && [ "$CMD" != "bash" ]; then
        dep1="-d afterok:${JOBIDS[${#JOBIDS[@]}-2]}"
        dep2="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    fi

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_1_fastqc.err \
        -o $dirList[$step]/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.raw1 \
        $dep1 runFastQC.sh $dirList[$step] $fastq1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_2_fastqc.err \
        -o $dirList[$step]/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.raw2 \
        $dep2 runFastQC.sh $dirList[$step] $fastq2`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

## Run Sortmerna. Depends on a successful run of FastQValidator
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    ## we do not depend on the FASTQC, but on FastQValidator if it
    ## was started
    if [ $pstart -eq (($step-1)) ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-3]}"
    fi
    JOBCMDS+=(`prepare_$CMD \
        $debug_var \
        -e $dirList[$step]/${sname}_sortmerna.err \
        -o $dirList[$step]/${sname}_sortmerna.out \
        $dep \
        -J ${sname}.RNAseq.SortMeRNA \
        runSortmerna.sh $sorterna_inx $dirList[$step] $fastq1 $fastq2`)
    if [ "$CMD" == "bash" ]; then
	JOBCMDS+=("export SORTMERNADIR=$SORTMERNADIR;$sortmerna_id")
    else
	JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
    fi

fi

fastq_sort_1="$dirList[$step]/${sname}_sortmerna_1.fq.gz"
fastq_sort_2="$dirList[$step]/${sname}_sortmerna_2.fq.gz"

# Run FastQC. Depends on a successful run of SortMeRNA
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -lt (($step+1)) ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    elif ([ ! -f $fastq_sort_1 ] || [ ! -f $fastq_sort_2 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_1_fastqc.err \
        -o $dirList[$step]/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.SortMeRNA1 \
        $dep runFastQC.sh $dirList[$step] $fastq_sort_1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)

    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_2_fastqc.err \
        -o $dirList[$step]/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.SortMeRNA2 \
        $dep runFastQC.sh $dirList[$step] $fastq_sort_2`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

## Run trimmomatic. Depends on a successful run of SortMeRNA
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -le (($step-1))  ] && [ "$CMD" != "bash" ]; then
        # SortMeRNA has to finish, if it was started
	# hence we check for a pstart lower or equal to 3
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-3]}"
    elif ([ ! -f $fastq_sort_1 ] || [ ! -f $fastq_sort_2 ]) && [ "$CMD" != "bash" ]; then
        # If we start here, make sure the files exist
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    JOBCMDS+=($(prepare_$CMD \
        $debug_var \
        -e $dirList[$step]/${sname}_trimmomatic.err \
        -o $dirList[$step]/${sname}_trimmomatic.log \
        -J ${sname}.RNAseq.Trimmomatic \
        $dep \
        runTrimmomatic.sh $phred_value $fastq_sort_1 $fastq_sort_2 $dirList[$step] $trimmomatic_options))
    JOBIDS+=($(run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}))
fi

## Trimmed fastq file paths
fastq_trimmed_1="$dirList[$step]/${sname}_sortmerna_trimmomatic_1.fq.gz"
fastq_trimmed_2="$dirList[$step]/${sname}_sortmerna_trimmomatic_2.fq.gz"

# Run FastQC. Depends on a successful run of Trimmomatic
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -lt (($step+1)) ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    elif ([ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_1_fastqc.err \
        -o $dirList[$step]/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.Trimmomatic1 \
        $dep runFastQC.sh $dirList[$step] $fastq_trimmed_1`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
    
    JOBCMDS+=(`prepare_$CMD \
        -e $dirList[$step]/${sname}_2_fastqc.err \
        -o $dirList[$step]/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.Trimmomatic2 \
        $dep runFastQC.sh $dirList[$step] $fastq_trimmed_2`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

# Run Salmon. Depends on a successful run of trimmomatic
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -le 5 ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-5]}"
    elif ([ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    
    JOBCMDS+=($(prepare_$CMD \
        $debug_var \
        -e $dirList[$step]/${sname}_salmon.err \
        -o $dirList[$step]/${sname}_salmon.out \
        -J ${sname}.RNAseq.kallisto \
        $dep runSalmon.sh $salmon_index $fastq_trimmed_1 $fastq_trimmed_2 $dirList[$step]))
    JOBIDS+=($(run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}))
fi

# Run Kallisto. Depends on a successful run of trimmomatic
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -le (($step-2)) ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-5]}"
    elif ([ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
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
  
    JOBCMDS+=($(prepare_$CMD \
        $debug_var \
        -e $dirList[$step]/${sname}_kallisto.err \
        -o $dirList[$step]/${sname}_kallisto.out \
        -J ${sname}.RNAseq.kallisto \
        $dep runKallisto.sh $kallisto_strand_option $fastq_trimmed_1 $fastq_trimmed_2 $kallisto_index $kallisto_fasta $dirList[$step]))
    JOBIDS+=($(run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}))
fi

# Run STAR. Depends on a successful run of Trimmomatic.
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    ## Trimmomatic needs finishing if it was started
    ## hence, we check if step 5 was started  
    if [ $pstart -le (($step-3)) ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-3]}"
    elif ([ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]) && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi

    if [ ! -z $star_gtf ]; then
	    star_runner_options="-g $star_gtf"
    fi
    
    if [ ! -z $star_intron_max ]; then
      star_runner_options="$star_runner_options -m $star_intron_max"
    fi

    JOBCMDS+=(`prepare_$CMD \
            $debug_var \
            -e $dirList[$step]/${sname}_STAR.err \
            -o $dirList[$step]/${sname}_STAR.out \
            -J ${sname}.RNAseq.STAR \
            ${mem:+"-m" "$mem"} \
            $dep runSTAR.sh -l "$bam_memory" $star_runner_options $dirList[$step] $star_ref $genome_fasta $fastq_trimmed_1 $fastq_trimmed_2 $star_options`)
    JOBIDS+=(`run_$CMD ${JOBCMDS[${#JOBCMDS[@]}-1]}`)
fi

# Run HTSeq. Depends on a successful run of STAR
(($step+=1))
if [ $pstart -le (($step+1)) ] && [ $pend -ge (($step+1)) ]; then
    echo "# Preparing step (($step+1))"

    dep=
    if [ $pstart -lt $step ] && [ "$CMD" != "bash" ]; then
        dep="-d afterok:${JOBIDS[${#JOBIDS[@]}-1]}"
    elif [ ! -f $star/${sname}_sortmerna_trimmomatic_STAR.bam ] && [ "$CMD" != "bash" ]; then
        echo >&2 "ERROR: STAR alignment BAM file could not be found"
        cleanup
    fi

    strand_arg=
    if [ $ilm_stranded -eq 1 ]; then
        strand_arg="-s"
    fi
    if [ $non_ilm_stranded -eq 1 ]; then
        strand_arg="-s -a"
    fi
    
    # HTSeq
    JOBCMDS+=(`prepare_$CMD \
        $debug_var \
        -e $dirList[$step]/${sname}_HTSeq.err \
        -o $dirList[$step]/${sname}_HTSeq.out \
        -J ${sname}.RNAseq.HTSeq \
        $dep runHTSeq.sh -i $idattr $strand_arg $dirList[$step] $star/${sname}_sortmerna_trimmomatic_STAR.bam $htseq_gff`)
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
