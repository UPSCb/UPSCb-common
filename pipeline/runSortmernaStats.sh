#!/bin/bash

usage(){
    echo >&2 \
"
    runSortmernaStats.sh <sortmerna dir>

    Arguments:
        sortmerna dir - the directory containing the SortMeRna reports

    Note:
        The UPSCb Environment Variable needs to be set to your
        Git UPSCb checkout dir.

    Details:
        It reads the SortMeRna log files and create a text file containing a 
        table to be readily added to the wiki.
"
exit 1
}

if [ $# -ne 1 ]; then
    echo "This function takes one argument"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the directory containing the SortMeRna reports"
    usage
fi

if [ -z $UPSCb ]; then
    echo "The UPSCb environment variable needs to be set."
    usage
fi

if [ ! -f $UPSCb/pipeline/runSortmernaStats.sh ]; then
    echo "Either your UPSC env. var. is not set correctly or your checkout is too old."
    usage
fi

if ! grep -i 'mtssu' $1/*.log >/dev/null 2>&1; then
    echo >$1/sortmernaStats.txt \
"|Sample|rRNA filtered reads|rRNA filtered reads (percent)|total rRNA|5S|5.8S|Arc16S|Bac16S|18S|Arc23S|Bac23S|28S|
|----|----|----|----|----|----|----|----|----|----|"
else
    echo >$1/sortmernaStats.txt \
"|Sample|total rRNA|5S|5.8S|16S|18S|23S|28S|mtSSU|
|----|----|----|----|----|----|----|----|----|"
fi
#grep "%" $1/*.log | awk '{if($2=="%"){if(NR > 1){printf "\n"};smpl=$1;sub(/^.*\//, "", smpl);sub("_sortmerna.log:","",smpl);printf "|"smpl"|"$4"%|"} else {printf $2"%|"}}' >> $1/sortmernaStats.txt

for f in $(find $1 -name "*.log"); do   
  echo "|"$(basename ${f/_sortmerna_rRNA.log//})"|"$(grep "%" $f | awk '{if($1=="Total"){if($3=="failing"){c2=$7;c3=$8;gsub(/\(|\)/,"",c3)}else{c4=$8;gsub(/\(|\)/,"",c4)}}else{c5t12=c5t12"|"$2}}END{print c2"|"c3"|"c4 c5t12"|"}') >> $1/sortmernaStats.txt
done
