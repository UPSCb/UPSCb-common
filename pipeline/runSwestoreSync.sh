#!/bin/bash -l
#SBATCH -p core -n 1
#SBATCH -A snic2018-13-9
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nicolas.delhomme@umu.se
#SBATCH -t 3-00:00:00

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# failsafe
set -ex

USAGETXT=\
"
	Usage: runSwestoreSync.sh <parent dir> <archive dir> <swestore project>
	
	Note: <archive dir> is relative to the parent dir!
"

# check arguments
if [ $# -ne 3 ]; then
    echo "This function needs 3 argument"
    usage
fi

if [ ! -d $1 ]; then
  abort "The first argument needs to be an existing directory"
fi

if [ ! -d $1/$2 ]; then
   abort "The second argument needs to be an existing directory within the first argument directory"
fi

# create a proxy
arcproxy -c validityPeriod=96H -p key=file:~delhomme/.globus/key.txt

# create the directory in swestore
dir=$(echo $2 | sed "s:^./::")
#arcmkdir $3/$dir

# go to the parent dir
cd $1

# create an archive - TODO check for special chars - checked manually, there are none
arx=$(echo $2 | sed 's:^./::' | sed 's:-:_:gm' | sed 's:/:-:').tar
tar -cf $arx $1/$dir

# sync the archive
arccp $arx $3/$dir/

# check
lsize=$(ls -l $arx | awk '{print $5}')
rsize=$(arcls -l $3/$dir/$arx | awk '{if(NR>1)print $3}')

if [ $lsize -ne $rsize ]; then
  abort "The file do not have the same size!"
else
  # remove the local structure and the archive
  rm $arx
  rm -rf $dir
fi
