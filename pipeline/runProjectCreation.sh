#!/bin/bash
set -e

# parameters
PI=
GID=
DESCRIPTION=
read -r -a KNOWNPIS <<< $(getent group | awk -F: '{if($3 >= 2007 && $3 < 3000){print $1}}' | grep -v u20 | xargs )
read -r -a KNOWNUIDS <<< $(getent passwd | awk -F: '{if($3 >= 1000 && $3 < 21000){print $1}}' | sort | xargs )
read -r -a KNOWNGROUPS <<< $(getent group | awk -F: '{if($3 >= 2007 && $3 < 3000){print $1}}' | grep u20 | xargs)

if [ -z $UPSCb ]; then
  export UPSCb=/mnt/picea/home/delhomme/Git/UPSCb
fi

# usage text
export USAGETXT=\
"
Usage: $0 <PI> <GID> <DESCRIPTION> <USER> [USER] [USER] ...

Details:

<PI> is one of:
${KNOWNPIS[@]}

<GID> is the new project ID: u20YYXXX with YY the year and
XXX the serial in that year
    
<DESCRIPTION> is the project description, make sure to 
quote it if you use white spaces. Avoid special characters
  
<USER> at least one user (known to LDAP) is mandatory

[USER] additional users, as many as required

Specifics:
  This script needs to be run on microasp
  You need to have the sudo privilege
  You need to have the LDAP admin password handy
"

# load functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../src/bash/functions.sh

# check the host
if [ $HOSTNAME != "microasp" ]; then
  abort "Not running on the right host"
fi

# get the arguments
if [ $# -lt 4 ]; then
  abort "The script expects a minimum of 3 arguments."
fi

# check the arguments
PI=$1
shift
if [ $(containsElement $PI "${KNOWNPIS[@]}") -eq 1 ]; then
  abort "Unknown PI"
fi

GID=$1
shift
if [ $(containsElement $GID "${KNOWNGROUPS[@]}") -eq 0 ]; then
  abort "This project ID already exists"
fi

DESCRIPTION=$1
shift

for u in "${@}"; do
  if [ $(containsElement $u "${KNOWNUIDS[@]}") -eq 1 ]; then
  abort "The user $u does not exist"
  fi
done

## LDAP
# create the group
sudo ldap_addgroup.sh $GID $DESCRIPTION

echo "Created the group: $GID"

# add the users
for u in "${@}"; do
  sudo ldap_addusertogroup.sh $GID $u
  echo "Added the user $u to the $GID group"
done

## SLURM
# create the account
sudo sacctmgr add account Description=$DESCRIPTION Organization=$PI  Name=$GID Parent=$PI
echo "Created the account $GID for the PI $PI"

# add the users
for u in "${@}"; do
  sudo sacctmgr add user Name=$u Accounts=$GID
  echo "Added the user $u to the $GID account"
done

# done
exit 0
