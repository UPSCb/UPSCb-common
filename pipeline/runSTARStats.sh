#!/bin/bash

usage(){
    echo >&2 \
"
    runSTARStats.sh <STAR dir>

    Arguments:
        STAR dir - the directory containing the STAR logs directory

    Note:
        The UPSCb Environment Variable needs to be set to your
        Git UPSCb checkout dir.

    Details:
        It reads the STAR 'Log.Final.out' log files and create a text file containing a 
        table to be readily added to the wiki.
"
exit 1
}

if [ $# -ne 1 ]; then
    echo "This function takes one argument"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the directory containing the STAR reports"
    usage
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    usage
fi

if [ ! -f $UPSCb/pipeline/runSTARStats.sh ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    usage
fi

echo > $1/STARStats.txt \
"|Sample|uniquely mapped|mismatch rate|deletion rate|insert rate|multiple mapping|too many mapping|unmapped MM|unmapped short|unmapped other|
|----------|---|---|---|---|---|---|---|---|"
cd $1
grep "%" */*Log.final.out | awk 'BEGIN{FS="|"}{pct=$2;gsub(" |\t","",pct);if (index($1,"Uniquely")>0){smpl=$1;sub("Log.final.out.*","",smpl);sub(".*/","",smpl);if (NR>1) {printf "\n"} printf"|"smpl"|"pct"%|"}else{printf pct"%|"}}' >> STARStats.txt

