#!/bin/bash

set -e

s1=$1
in=$2

#Argument number and type tests
if [ $# != 2 ]; then
	echo "This script takes 2 arguments, an identifier for samples that are to be merged
	and an input directory"
	exit 1
fi

if [ ! -d $in ]; then
	echo "Your input argument is not a directory"
	exit 1
fi
#Find relevant samples and decompress
find $in -name "*$s1*" | grep sortmerna | grep $s1 | xargs -I {} -P 3 gzip -d {}
#Skip (only rename) if files to be merged are less than 2, merge others
#The process is gzip -d -> find reads -> merge (with cat) -> gzip
number=`find $in -name "*.fq" | grep $s1 | grep "_1" | wc -l`
if [ $number -lt 2 ];then
	echo "SKIPPING but renaming"
	fw=`find $in -name "*.fq" | grep $s1 | grep "_1" | head -1`
	echo "$fw"
	mv $f ${f/_1.fq/}_1_merged.fq
	echo "AND"
	rv=`find $in -name "*.fq" | grep $s1 | grep "_2" | head -1`
        echo "$rv"
	mv $f ${f/_1.fq/}_1_merged.fq

else
	echo "MERGING FORWARD"
	find $in -name "*.fq" | grep $s1 | grep "_1"
	echo "TO:"
	f=`find $in -name "*.fq" | grep $s1 | grep "_1" | head -1`
	echo ${f/_1.fq/}_1_merged.fq
	cat `find $in -name "*.fq" | grep $s1 | grep "_1"` > ${f/_1.fq/}_1_merged.fq
	gzip ${f/_1.fq/}_1_merged.fq
	#place cat here
	echo "AND REVERSE"
	find $in -name "*.fq" | grep $s1 | grep "_2"
	echo "TO:"
        echo ${f/_1.fq/}_2_merged.fq
	cat `find $in -name "*.fq" | grep $s1 | grep "_1"` > ${f/_1.fq/}_2_merged.fq
        gzip ${f/_1.fq/}_2_merged.fq
	#pl
fi
#Clean up behind
find $in -name "*$s1*" | grep sortmerna | grep $s1 | xargs -I {} -P 3 gzip {}
