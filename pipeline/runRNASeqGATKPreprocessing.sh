#!/bin/bash -l

set -ex

if [ $# -lt 5 ]; then
    echo >&2 "
      Usage: $0 [options] <proj> <mail> <bam> <ref.fasta> <output directory>
    
    Options:
    
      -d dbSNP file; implies that BaseRecalibration is run
      -f forces the whole pipeline to be rerun
      -m define the minimal mate distance 
      -p turns GATK --fix_misencoded_quality_scores to TRUE  
    "
    exit 1
fi

# OPTION BLOCK
OPTIONS=
MIN=-1
FORCE=0
dbSNP=
while getopts "d:fm:p" opt; do
    case "$opt" in
        d) dbSNP=$OPTARG ;;
        f) FORCE=1 ;;
        m) MIN=$OPTARG ;;
        p) OPTIONS="--fix_misencoded_quality_scores" ;;
    esac
done
shift `expr $OPTIND - 1`

proj=$1
mail=$2
inbam=$3
ref=$4
outdir=$5

if [ ! -f $inbam ]; then
    echo >&2 "bam file not found: $inbam"
    exit 1
fi

if [ ! -f $ref ]; then
    echo >&2 "reference fasta file not found: $ref"
    exit 1
fi

if [ ! -z $dbSNP ] && [ ! -f $dbSNP ]; then
    echo >&2 "dbSNP file not found: $dbSNP"
    exit 1
fi

# Sample name
bname=$(basename $inbam)
#sname=${bname%%_*}
sname=${bname//.bam/}

# Create directory structure
[[ ! -d $outdir ]] && mkdir $outdir

markdupdir="$outdir/markduplicates"
[[ ! -d $markdupdir ]] && mkdir $markdupdir

splittrimdir="$outdir/splitntrim"
[[ ! -d $splittrimdir ]] && mkdir $splittrimdir

if [ ! -z $dbSNP ]; then
  basrecaldir="$outdir/baserecalibration"
  [[ ! -d $basrecaldir ]] && mkdir $basrecaldir
fi

indelrealigndir="$outdir/indelrealign"
[[ ! -d $indelrealigndir ]] && mkdir $indelrealigndir

PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"

run_sbatch_usage() {
    echo "usage: run_sbatch OPTIONS <batch script> [<script param> ...]

    Run a batch script and echo the job id

    OPTIONS
      -e  stderr (required)
      -o  stdout (required)
      -J  job name
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

    while getopts "J:e:o:d:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            ?) run_sbatch_usage; exit 1 ;;
        esac
    done

    shift $((OPTIND-1))

    script=$1
    shift

    if [ -z $jobname ]; then
        jobname="${sname}.${script}"
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

    sbatch_echo=`sbatch -A "$proj" \
                -J "$jobname" \
                --mail-user "$mail" \
                --mail-type=ALL \
                -e "$log_path" \
                -o "$out_path" \
                $sbatch_options \
                $PIPELINE_DIR/$script $@`

    if [ $? -ne 0 ]; then
        echo "ERROR: submission failed" 1>&2
        cleanup ${JOBIDS[*]}
    fi

    echo ${sbatch_echo//[^0-9]/}
}

JOBIDS=()

# Run MarkDuplicates
markdup_bam=$markdupdir/$(basename ${inbam/.bam/})_mkdup.bam
if [ ! -f $markdup_bam ]; then
    markdup_id=$(run_sbatch \
        -e $markdupdir/${sname}_markdup.err \
        -o $markdupdir/${sname}_markdup.out \
        runMarkDuplicates.sh $inbam $markdupdir)
    JOBIDS+=($markdup_id)
fi

# Run SplitNCigarReads
splittrim_bam=$splittrimdir/$(basename ${markdup_bam/.bam/})_split.bam
if [ ! -f $splittrim_bam ]; then
    dep=
    [[ ! -f $markdup_bam ]] && dep="-d afterok:$markdup_id"
    splittrim_id=$(run_sbatch \
        -e $splittrimdir/${sname}_split.err \
        -o $splittrimdir/${sname}_split.out \
        $dep \
        runGATK_SplitNCigarReads.sh $markdup_bam $ref $splittrimdir)
    JOBIDS+=($splittrim_id)
fi

# TODO add a block here for bqsr
# Run BaseRecalibrator
if [ ! -z $dbSNP ]; then

  bqsr_targets=$basrecaldir/$(basename ${splittrim_bam/.bam/_recalibrated.bam})
  
  if [ ! -f $bqsr_targets ] || [ $FORCE -eq 1 ]; then
    dep=
    [[ ! -f $splittrim_bam ]] && dep="-d afterok:$splittrim_id"
    bqsr_target_id=$(run_sbatch \
	    -e $basrecaldir/${sname}.err \
      -o $basrecaldir/${sname}.out \
	    $dep \
	    runGatk_BaseRecalibration.sh $splittrim_bam $ref $basrecaldir $dbSNP)
    JOBIDS+=($bqsr_target_id)
  fi
else
  bqsr_targets=$splittrim_bam
  bqsr_target_id=$splittrim_id
fi


# Run RealignerTargetCreator
indelrealign_intervals=$indelrealigndir/$(basename ${bqsr_targets/.bam/.intervals})
if [ ! -f $indelrealign_intervals ]; then
    dep=
    [[ ! -f $bqsr_targets ]] && dep="-d afterok:$bqsr_target_id"
    indelinterval_id=$(run_sbatch \
        -e $indelrealigndir/${sname}_interval.err \
        -o $indelrealigndir/${sname}_interval.out \
        $dep \
        runGatkRealignerTargetCreator.sh $bqsr_targets $ref $indelrealigndir)
    JOBIDS+=($indelinterval_id)
fi

# Run IndelRealigner
indelrealign_bam=$indelrealigndir/$(basename ${indelrealign_intervals/.intervals/.bam})
if [ ! -f $indelrealign_bam ]; then
    dep=
    [[ ! -f $bqsr_targets ]] && dep="-d afterok:$indelinterval_id"
    indelrealign_id=$(run_sbatch \
        -e $indelrealigndir/${sname}_realign.err \
        -o $indelrealigndir/${sname}_realign.out \
        $dep \
        runGATK_IndelRealigner.sh $bqsr_targets $ref $indelrealign_intervals $indelrealigndir)
    JOBIDS+=($indelrealign_id)
fi

echo "Successfully submitted ${#JOBIDS[@]} jobs for ${sname}: ${JOBIDS[@]}"
