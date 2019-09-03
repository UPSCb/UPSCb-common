#!/bin/bash -l
#SBATCH --mem=100G
#SBATCH --mail-type=all
#SBATCH -n 1
#SBATCH -J scThresh

#Make sure we have igraph
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

usage(){
    echo "Usage: runThresholdNetwork.sh <input file> <output file> ['other options']"
    echo "Important: Other options must be quoted"

    exit 1
}

INF=$1
OUTF=$2
OPTS=""

shift 2
if [ $# == 1 ]; then
    OPTS=$1
fi

if [ $# -gt 1 ];then
    usage
fi

if [ ! -f $INF ];then
    usage
fi

if [ -z $UPSCb ];then
    echo "You must set the UPSCb environment variable for this script"
    usage
fi

$UPSCb/src/cpp/scgraph/threshold -i $INF $OPTS > $OUTF
