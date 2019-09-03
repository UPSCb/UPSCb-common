#!/bin/bash -l
#SBATCH -A b2010064
#SBATCH -n 8
#SBATCH -t 1-0:00:00
#SBATCH --mail-user david.sundell@plantphys.umu.se
#SBATCH --mail-type=ALL

prog="/home/davidsu/opt/Trimmomatic-0.22/trimmomatic_0.22.sh"

####
##		Run trimmomatic
###

## Usage: runTrim.sh inptFolder forwardReads reverseReads

cd $1
OUT="/proj/b2010064/nobackup/david_pipeline/$4/trimmed"
n1=${2%.fq*}
n2=${3%.fq*}
in1=$2
in2=$3
if [ -h $f ]; then
	in1=$(readlink -f "$2")
	in2=$(readlink -f "$3")
fi

#run
sh $prog $in1 $in2 $OUT/$n1"_FP.fq.gz" $OUT/$n1"_FU.fq.gz" $OUT/$n2"_FP.fq.gz" $OUT/$n2"_FU.fq.gz" ILLUMINACLIP:"/proj/b2010064/analysis/illuminaClippingPoly.fa":2:40:14 MINLEN:50