#!/bin/bash
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 20
## time too for large files
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

## stop on error
set -e

## be verbose and extend the commands
set -x

## check the options if any
#KEEP=1
#useMtSSU=0
UNPAIRED=0
PROC=20
DBS=
SEEDS=2
PASSES=
SAM="--sam --num_alignments 1"
ILV=0

## source functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

## usage
export USAGETXT="
	Usage: runSortmerna.sh [option] <out dir> <tmp dir> <forward fastq.gz> <reverse fastq.gz>

	Options:
                -a report alignments (default is on, set to skip)
                -d define your dbs (semi-colon separated)
#                -k drop the rRNA (only for v1.9, default to keep them)
#                -m run against mtSSU in addition (only for v1.9)
                -p number of threads to be used (default $PROC)
                -P number of passes (default to sortmerna defaults)
		            -s number of seeds (default 2)
                -u single end data (in that case only the forward fastq is needed)
		-i interleaved data (in that case only one fastq is needed)

         Note:
               1) The SORTMERNA_DBS environment variable needs to be set
               2) Only SortMeRna version 2 and higher are supported
#               3) -m is not applicable if -d is set
               4) the 2nd argument is just a place holder on UPPMAX, it is automatically replaced by the node local tmp
"

## then check for availability
# tool=`which sortmerna 2>/dev/null`
# if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
#   echo "sortmerna available"
# else
#   echo "ERROR: INSTALL SortMeRna"
#   usage
# fi

## check for sortmerna version
version=$(sortmerna --version 2>&1 | grep "version" | cut -d" " -f3 | cut -d. -f 1)
#is1dot9=`sortmerna --version 2>&1 | grep version | grep 1.9 | wc -c`
#is2dotx=`sortmerna --version 2>&1 | grep "version 2." | wc -c`

#if [ $is1dot9 == 0 ] && [ $is2dotx  == 0 ]; then
if [ $version -lt 2 ]; then
  abort "Only version 2 and higher are supported"
fi

## get the options
#while getopts d:kmp:u option
while getopts ad:ip:P:s:u option
do
        case "$option" in
        a) SAM=;;
      d) DBS=$OPTARG;;
	i) ILV=1
	   UNPAIRED=1;;
#	    k) KEEP=0;;
#	    m) useMtSSU=1;;
	    p) PROC=$OPTARG;;
	    P) PASSES="--passes $OPTARG";;
	    s) SEEDS=$OPTARG;;
	    u) UNPAIRED=1;;
		\?) ## unknown flag
		abort;;
        esac
done
shift `expr $OPTIND - 1`

##
echo Setting up

## set some env var
## this location is not in Git anymore!
## it has to be downloaded by the user
## check the ethylene-insensitive project submitter to see
## how to set that up
if [ -z $SORTMERNA_DBS ]; then
    abort "You need to set your SORTMERNA_DBS environment variable"
fi

## set the default dbs
if [ ! -z $DBS ]; then
  dbs=${DBS//;/ }
  dbNum=`echo $DBS | awk -F";" '{print NF}'`
else
#  if [ $is2dotx != 0 ]; then
    db5s=$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98
    db58s=$SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98
    db16sa=$SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95
    db16s=$SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90
    db18s=$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95
    db23sa=$SORTMERNA_DBS/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DBS/index/silva-arc-23s-id98
    db23s=$SORTMERNA_DBS/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DBS/index/silva-bac-23s-id98
    db28s=$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98
    dbs="$db5s:$db58s:$db16sa:$db16s:$db18s:$db23sa:$db23s:$db28s"
#  else
#    db5s=$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta
#    db58s=$SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta
#    db16sa=$SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta
#    db16s=$SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id85.fasta
#    db18s=$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta
#    db23sa=$SORTMERNA_DBS/rRNA_databases/silva-arc-23s-id98.fasta
#    db23s=$SORTMERNA_DBS/rRNA_databases/silva-bac-23s-id98.fasta
#    db28s=$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta
#    dbNum=8
#    dbs="$db5s $db58s $db16sa $db16s $db18s $db23sa $db23s $db28s"
#  fi

#  ## Add the mtSSU
#  if [ $is1dot9 != 0 ] && [ $useMtSSU == 1 ]; then
#    mtSSU=$SORTMERNA_DBS/rRNA_databases/mtSSU_UCLUST-95-identity.fasta
#    dbs="$dbs $mtSSU"
#    dbNum=9
#  fi
fi

##
echo Checking

## we get two dir and two files as input
if [ $UNPAIRED == 0 ]; then
    if [ $# != 4 ]; then
	abort "This function takes two directories and two files as arguments"
    fi
else
    if [ $# != 3 ]; then
	abort "This function takes two directories and one file as argument"
    fi
fi

if [ ! -d $1 ]; then
    abort "The first argument needs to be an existing directory"
fi

if [ ! -d $2 ]; then
    abort "The second argument needs to be an existing directory"
fi

## UPPMAX hack
if [ ! -z $SNIC_RESOURCE ]; then
  tmp=/scratch/$SLURM_JOB_ID
  # even uglier hack; should be fixed as soon as uppmax add the merge scripts to the module file
  export PATH=/sw/apps/bioinfo/SortMeRNA/2.1b/src_milou/sortmerna-2.1b/scripts:$PATH
else
  tmp=$2
fi
echo TMP: $tmp

## 
echo Gunzipping

## unzip the files
if [ ! -f $3 ]; then
    abort "The third argument needs to be an existing fastq.gz file"
fi
f1=`basename ${3//.gz/}`

if [ $UNPAIRED == 0 ]; then
    if [ ! -f $4 ]; then
	abort "The forth argument needs to be an existing fastq.gz file"
    fi
    f2=`basename ${4//.gz/}`
fi

## decompress them
if [ ! -f $tmp/$f1 ]; then
    gunzip -c $3 > $tmp/$f1
fi
if [ $UNPAIRED == 0 ]; then
    if [ ! -f $tmp/$f2 ]; then
	gunzip -c $4 > $tmp/$f2
    fi
fi

## interleave them
fm=`basename ${3//.f*q.gz/}`
if [ $UNPAIRED == 0 ]; then
  merge-paired-reads.sh $tmp/$f1 $tmp/$f2 $tmp/$fm
fi

##
if [ $UNPAIRED == 0 ]; then
    echo Pre-cleaning
    rm -f $tmp/$f1 $tmp/$f2
else
    echo "TODO: Cleaning needs implementing for single end sequencing"
fi

##
echo Sorting

## PE
if [ $UNPAIRED == 0 ]; then
    fo=`basename ${3//_[1,2].f*q.gz/_sortmerna}`
else
    fo=`basename ${3//.f*q.gz/_sortmerna}`
fi

## check the options
opt="-a $PROC"

#if [ $KEEP == 1 ] && [ $is1dot9 != 0 ]; then
#  opt="$opt --bydbs --accept $tmp/${fo}_rRNA"
#fi 

## run
if [ $UNPAIRED == 0 ]; then
#  if [ $is2dotx != 0 ]; then
    sortmerna --ref $dbs --reads $tmp/$fm --other $tmp/$fo --log --paired_in --fastx $opt \
$SAM $PASSES --num_seeds $SEEDS --aligned $tmp/${fo}_rRNA
#  else
#    sortmerna -n $dbNum --db $dbs --I $tmp/$fm --other $tmp/$fo --log $1/$fo --paired-in $opt
#  fi
else
#  if [ $is2dotx != 0 ]; then
    if [ $ILV -eq 1 ]; then
       sortmerna --ref $dbs --reads $tmp/$f1 --other $tmp/$fo --log --paired_in --fastx $opt \
$SAM $PASSES --num_seeds $SEEDS --aligned $tmp/${fo}_rRNA
       UNPAIRED=0
    else 
       sortmerna --ref $dbs --reads $tmp/$f1 --other $1/$fo --log $opt $SAM \
--fastx $PASSES --num_seeds $SEEDS --aligned $tmp/${fo}_rRNA
    fi
#  else
#    sortmerna -n $dbNum --db $dbs --I $tmp/$f1 --other $1/$fo --log $1/$fo $opt
#  fi
fi

## deinterleave it
if [ $UNPAIRED == 0 ]; then
    ## sortmerna get confused by dots in the filenames
    if [ ! -f $tmp/$fo.fastq ]; then
	    mv $tmp/$fo.* $tmp/$fo.fastq
    fi
    unmerge-paired-reads.sh $tmp/$fo.fastq $1/${fo}_1.fq $1/${fo}_2.fq
fi

## cleanup
echo Post-Cleaning

#if [ $is2dotx != 0 ]; then
  ## mv the rRNA, fastq and log back
  mv $tmp/${fo}_rRNA.* $1
#fi

## rm the tmp
if [ $UNPAIRED == 0 ]; then
    rm -f $tmp/$fo.fastq
    if [ $ILV -eq 1 ]; then
      rm -f $tmp/$f1
    else
      rm -f $tmp/$fo
    fi 
else
    rm -f $tmp/$f1
fi

## deinterleave the rest if needed
#if [ $KEEP == 1 ]; then
#    if [ $UNPAIRED == 0 ]; then
#	find $tmp -name "${fo}_rRNA*" -print0 | xargs -0 -I {} -P 6 sh -c 'unmerge-paired-reads.sh $0 $1/`basename ${0//.fastq/_1.fq}` $1/`basename ${0//.fastq/_2.fq}`' {} $1
#    fi
#fi

## keep that as a reminder if that happens again
## sortmerna get confused by the dots as well...
## echo Validating
if [ $UNPAIRED == 1 ]; then
    if [ ! -f $1/$fo.fastq ] && [ ! -f $1/$fo.fq ] ; then
        abort "Could not find the output file. Check the files: $1/$fo"
        #find $1 -name "$fo*" | grep -v *.err | grep -v *.out | xargs
        #mv $1/$fo.$1/$fo.fq
    fi
fi

## 
echo Gzipping

## compress the output files
if [ $UNPAIRED == 0 ]; then
  find $1 -name "${fo}*.fq" -print0 | xargs -0 -I {} -P $PROC gzip -f {}
else
  if [ -f $1/$fo.fastq ]; then
    gzip -c $1/$fo.fastq > $1/$fo.fq.gz
    rm $1/$fo.fastq
  else
    gzip $1/$fo.fq
  fi
fi

# compress the 2.X new files (sam and _rRNA)
#if [ $is2dotx != 0 ]; then
  find $1 -name "${fo}_rRNA.[f,s]*" -print0 | xargs -0 -I {} -P $PROC gzip -f {}
#fi

##
echo Done
