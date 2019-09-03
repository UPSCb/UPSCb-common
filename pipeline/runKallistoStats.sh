#!/bin/bash

usage(){
    echo >&2 \
"
    runkallistoStats.sh <kallisto dir>

    Arguments:
        kallisto dir - the directory containing the kallisto stderr files

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
    echo "The first argument needs to be the directory containing the kallisto reports"
    usage
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    usage
fi

if [ ! -f $UPSCb/pipeline/runKallistoStats.sh ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    usage
fi

cd $1
grep pseudoaligned *.err | awk '{smpl=$1;gsub(/_kallisto.*$/,"",smpl);val=$5;gsub(/,/,"",val);print smpl,val}' > kallistoStats.txt
