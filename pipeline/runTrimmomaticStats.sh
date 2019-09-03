#!/bin/bash

usage(){
    echo >&2 \
"
    runTrimmomaticStats.sh <trimmomatic dir>

    Arguments:
        trimmomatic dir - the directory containing the Trimmomatic error logs

    Note:
        The UPSCb Environment Variable needs to be set to your
        Git UPSCb checkout dir.

    Details:
        It reads the Trimmomatic error log files and create a text file containing a 
        table to be readily added to the wiki.
"
exit 1
}

if [ $# -ne 1 ]; then
    echo "This function takes one argument"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the directory containing the Trimmomatic reports"
    usage
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    usage
fi

if [ ! -f $UPSCb/pipeline/runTrimmomaticStats.sh ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    usage
fi

echo > $1/trimmomaticStats.txt \
"|Sample|Both|B%|Forward|F%|Reverse|R%|Dropped|D%|
|----|----|----|----|----|----|----|----|----|"
grep "Input" $1/*.err | awk -F" " '{gsub(/\(/,"");gsub(/\)/,"");smpl=$1;sub(/^.*\//, "", smpl);sub("_trimmomatic.err:Input","",smpl);print "|"smpl"|"$7"|"$8"|"$12"|"$13"|"$17"|"$18"|"$20"|"$21"|"}' >> $1/trimmomaticStats.txt

