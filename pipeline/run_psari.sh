#! /bin/bash -l

echo -n "Time started: " 
date

set -ex

#SBATCH -A b2016040
#SBATCH -o /mnt/picea/home/mdong/Git/UPSCb/projects/spruce-pseudogene/reports/psari/pseudo.out
#SBATCH -e /mnt/picea/home/mdong/Git/UPSCb/projects/spruce-pseudogene/reports/psari/pseudo.err
#SBATCH -J psari.job
#SBATCH -p core
#SBATCH -c 8
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mickael.dong@u-psud.fr
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load perl

TmpID=${SLURM_JOB_ID:-$$}

BlastInput=$1
RunDir=/mnt/picea/projects/spruce/pipeline/psari_data

Prefix=${BlastInput%.txt.gz}.psari

echo -n "Time started psari.pl: " 
date

DiscoveryOutput=$Prefix.prelim.gff3

#scripts_dir=/mnt/picea/Modules/apps/bioinfo/psari/16.04.13/bin
#cd /mnt/picea/Modules/apps/bioinfo/psari/16.04.13/bin
scripts_dir=/mnt/picea/home/mdong/Git/pseudogene/pipelines/psari

perl $scripts_dir/psari.pl $BlastInput > $DiscoveryOutput

echo -n "Time started psari_remove-HC-overlaps.sh: " 
date

ScreeningOutput=$Prefix.HC-screened.gff
ScreeningOverlapOutput=$Prefix.HC-overlap.gff  # those that do overlap HC genes

sh $scripts_dir/psari_screen-HC-overlaps.sh $DiscoveryOutput

# The alternate screening script is perhaps more complete but may not run because
# it includes a grep that may gobble up too much memory
# ./psari_screen-HC-overlaps_alternate.sh $DiscoveryOutput  

echo -n "Time started psari_annotateBedWithVariants: " 
date

PsariOutput=$Prefix.full.gff3

sh $scripts_dir/annotateBedWithVariants.sh $RunDir/$ScreeningOutput | gzip -c > $RunDir/$PsariOutput.gz

# convenience scripts produced by the annotate script
rm -f psari_join-annotation.pl psari_produce-annotation.awk  


echo -n "Time ended: " 
date
echo
echo "Output in $RunDir/$PsariOutput.gz" 

