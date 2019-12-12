#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:01:00
#SBATCH --mail-type=ALL

set -ex

usage(){
  echo >&2 \
  "
  Usage: $0 [options] <run dir>

  Note:
  This script expects one argument: the root directory where the RNA-Seq
  pipeline was run.

  Options:
    -t the directory to which the results should be transferred to using scp.
  "
  exit 1
}

# global var
SCP=

# get the options
while getopts t: option
do
  case "$option" in
	    t) SCP=$OPTARG;;
		  \?) ## unknown flag
		  echo "There is no such option"
		  usage;;
  esac
done
shift `expr $OPTIND - 1`

# check
if [ $# != 1 ]; then
  echo "ERROR: This script expects one argument"
  usage
fi

if [ ! -d $1 ]; then
  echo "ERROR: The first argument should be an existing directory"
  usage
fi

# setup
dir=$(realpath $1)
#dir=$1
cd $dir
mkdir -p reports
exec=$(realpath $(dirname $0))

# collate the fastqc report
# raw
cd $dir/fastqc/raw
$exec/runFastQCMultiviewer.sh .
echo -e "sample\trawReads" > $dir/reports/rawCounts.txt
find multiview/*_1_f* -name fastqc_data.txt | sort | xargs -I {} bash -c 'grep "Filename" $0 | awk "BEGIN{ORS=\"\t\"}{print \$2}"; grep "Total Sequences" $0 | awk "{print \$3}"' {} >> $dir/reports/rawCounts.txt
tar -zcf ../../reports/raw-multiview.tgz multiview multiview.html

# sortmerna
cd $dir/fastqc/sortmerna
$exec/runFastQCMultiviewer.sh .
echo -e "sample\tsortedReads" > $dir/reports/sortmernaCounts.txt
find multiview/*_1.fq* -name fastqc_data.txt | sort | xargs -I {} bash -c 'grep "Filename" $0 | awk "BEGIN{ORS=\"\t\"}{print \$2}"; grep "Total Sequences" $0 | awk "{print \$3}"' {} >> $dir/reports/sortmernaCounts.txt
tar -zcf ../../reports/sortmerna-multiview.tgz multiview multiview.html

# trimmomatic
cd $dir/fastqc/trimmomatic
$exec/pipeline/runFastQCMultiviewer.sh .
echo -e "sample\ttrimmedReads" > $dir/reports/trimmomaticCounts.txt
find multiview/*_1.fq* -name fastqc_data.txt | sort | xargs -I {} bash -c 'grep "Filename" $0 | awk "BEGIN{ORS=\"\t\"}{print \$2}"; grep "Total Sequences" $0 | awk "{print \$3}"' {} >> $dir/reports/trimmomaticCounts.txt
tar -zcf ../../reports/trimmomatic-multiview.tgz multiview multiview.html

# calculate the stats
# sortmerna
cd $dir/sortmerna
$exec/pipeline/runSortmernaStats.sh .
mv sortmernaStats.txt ../reports

# trimmomatic
cd $dir/trimmomatic
$exec/pipeline/runTrimmomaticStats.sh .
mv trimmomaticStats.txt ../reports

# STAR
cd $dir/star
echo -e "sample\taligned" > ../reports/star-counts.txt
find *_logs -name "*Log.final.out" | sort | xargs -I {} bash -c 'echo $0 | awk -F_ "{printf \"%s_%s\\t\",\$4,\$5}"; grep mapped $0 | grep -i number | grep -v many | awk "BEGIN{SUM=0}{SUM+=\$NF}END{print SUM}"' {} >> ../reports/star-counts.txt

$exec/pipeline/runSTARStats.sh .
mv STARStats.txt ../reports

# HTSeq
cd $dir/htseq
echo -e "sample\tnoFeature\tAmbiguous\tLowQual\tUnaligned\tnotUnique" > ../reports/htseq-counts.txt
find . -name "*.txt" | sort | xargs -I {} bash -c 'echo $0 | awk -F_ "{printf \"%s_%s\",\$4,\$5}"; tail -5 $0 | awk "{printf \"\\t%s\",\$2}"; echo' {} >> ../reports/htseq-counts.txt

# Kallisto
cd $dir/kallisto
$exec/pipeline/runKallistoStats.sh .
mv kallistoStats.txt ../reports

# create a final archive
cd $dir
tar -zcf reports.tgz reports

# export if the option is set
if [ ! -z $SCP ]; then
  scp reports.tgz $SCP
fi
