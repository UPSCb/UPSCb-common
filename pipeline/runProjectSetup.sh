#!/bin/bash
set -e

# parameters
PI=
SPECIES=
GID=
DATADIR=/mnt/picea/storage/data
PROJECTDIR=/mnt/picea/projects
read -r -a KNOWNPIS <<< $(getent group | awk -F: '{if($3 >= 2007 && $3 < 3000){print $1}}' | grep -v u20 | xargs )
read -r -a KNOWNSPECIES <<< $(find $DATADIR -mindepth 1 -maxdepth 1 -type d -exec basename "{}" \; | xargs )
read -r -a KNOWNGROUPS <<< $(getent group | awk -F: '{if($3 >= 2007 && $3 < 3000){print $1}}' | grep u20 | xargs )

if [ -z $UPSCb ]; then
  export UPSCb=/mnt/picea/home/delhomme/Git/UPSCb
fi

# usage text
export USAGETXT=\
"
Usage: $0 <SPECIES> <PI> <GID> <NAME>

Details:

<SPECIES> is one of:
${KNOWNSPECIES[@]}

<PI> is one of:
${KNOWNPIS[@]}

<GID> is the new project ID: u20YYXXX with YY the year and
XXX the serial in that year; one of:
${KNOWNGROUPS[@]}

<NAME> is the project name, without white space and special
characters

Specifics:
  This script needs to be run on picea
"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../src/bash/functions.sh

# check the host
if [ $HOSTNAME != "picea" ]; then
  abort "Not running on the right host"
fi

# get the arguments
if [ $# -ne 4 ]; then
  abort "The script expects 4 arguments."
fi

# check the arguments
SPECIES=$1
if [ $(containsElement $SPECIES "${KNOWNSPECIES[@]}") -eq 1 ]; then
  abort "Unknown species"
fi

PI=$2
if [ $(containsElement $PI "${KNOWNPIS[@]}") -eq 1 ]; then
  abort "Unknown PI"
fi

GID=$3
if [ $(containsElement $GID "${KNOWNGROUPS[@]}") -eq 1 ]; then
  abort "Unknown project ID"
fi

NAME=$4

# create the data directory
# and set up the permissions
DATADIR=$DATADIR/$SPECIES/$PI
if [ ! -d $DATADIR ]; then
  echo "Creating the directory: $DATADIR"
  sudo chattr -i $(dirname $DATADIR)
  sudo mkdir $DATADIR
  sudo chmod 771 $DATADIR
  sudo chgrp $GID $DATADIR
  sudo chmod g+s $DATADIR
  sudo chattr +i $(dirname $DATADIR)
fi

sudo chattr -i $DATADIR
sudo mkdir -p $DATADIR/$NAME
sudo chmod 771 $DATADIR/$NAME
sudo chgrp $GID $DATADIR/$NAME
sudo chmod g+s $DATADIR/$NAME
sudo chattr +i $DATADIR
sudo setfacl -Rm d:g:$GID:rwx,g:$GID:rwx $DATADIR/$NAME

echo "Created the directory: $DATADIR/$NAME"

# create the project directory
# and set up permissions
PROJECTDIR=$PROJECTDIR/$SPECIES/$PI
if [ ! -d $PROJECTDIR ]; then
  sudo chattr -i $(dirname $PROJECTDIR)
  sudo mkdir $PROJECTDIR
  sudo chattr +i $(dirname $PROJECTDIR)
fi

sudo mkdir -p $PROJECTDIR/$NAME/raw
sudo chmod -R 771 $PROJECTDIR/$NAME

sudo chgrp -R $GID $PROJECTDIR/$NAME
sudo chmod -R g+s $PROJECTDIR/$NAME
sudo setfacl -Rm d:g:$GID:rwx,g:$GID:rwx $PROJECTDIR/$NAME

echo "Created the directory: $PROJECTDIR/$NAME/raw"

# done
exit 0
