#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 4
#SBATCH -t 0-02:00:00
#SBATCH --mail-type=ALL

## stop on error and be verbose in the output
set -e -x



### tool sanity
if [ ! -z $SLURM_SUBMIT_DIR ]; then
    module load bioinfo-tools
    module load samtools/0.1.19
  
   
else
	samtools=`which samtools`
	if [ $? != 0 ]; then
		echo "please install samtools before running this script or add it to your PATH"
		exit 1
	fi

	if [ ! -f $samtools -a ! -x $samtools ]; then
		echo "your samtools does not appear to be an executable file"
		exit 1
	fi
fi

  



## arguments                                                                                                                                                                         
if [ $# != 2 ]; then
   echo "This script takes two arguments: the input file and the output directory"
   exit 1
fi

## input file                                                                                                                                                                        
if [ ! -f $1 ]; then
        echo "The first argument needs to be an existing file"
        exit 1
fi

## output dir                                                                                                                                                                        
if [ ! -d $2 ]; then
        echo "The second argument needs to be an existing output directory."
fi



## start           
samtools view $1 | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > $2/`basename $1`_counts.log


