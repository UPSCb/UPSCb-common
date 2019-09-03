#! /bin/bash
#SBATCH -p core -n 1
#SBATCH -t 24:00:00

# be verbose and stop on error
set -eux

# test 
#isEnvVarSet $UPSCb

# load module
module load R

if [ ! -d $1 ]; then
    echo "The first argument needs to be a directory"
    exit 1
fi
# exit 1 = exit with any error (normal case would be "exit 0" -> no error)


# one parameter from the command line: the saturation analysis directory. $1=1st argument on the command line

cd $1
Rscript --vanilla $UPSCb/projects/DAP-Seq/src/R/parseMacs2.R