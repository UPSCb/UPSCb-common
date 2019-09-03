#!/bin/bash

usage(){
    echo >&2 \
"
    runSalmonStats.sh <Salmon dir>

    Arguments:
        Salmon dir - the directory containing the Salmon stderr files

    Note:
        The UPSCb Environment Variable needs to be set to your
        Git UPSCb checkout dir.

    Details:
        It reads the kallisto stderr files and create a text delimited file
        containing the sample name and number of pseudoalignments
"
exit 1
}

if [ $# -ne 1 ]; then
    echo "This function takes one argument"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the directory containing the Salmon reports"
    usage
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    usage
fi

if [ ! -f $UPSCb/pipeline/runSalmonStats.sh ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    usage
fi

cd $1
grep Counted *.err | awk '{smpl=$1;gsub(/_sortmerna.*$/,"",smpl);val=$6;print smpl,val}' > SalmonStats.txt
