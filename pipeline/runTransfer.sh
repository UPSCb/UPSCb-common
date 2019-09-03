#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL

## stop on error
set -e

## helper function
errorState(){
    if [ $# -ge 1 ];then
	echo "ERROR: Exiting as $1"
	if [ ! -z $2 ];then
	        rm $2
		    checkStatus "the file removal failed. Do so manually!"
		    fi
    else
	echo "ERROR: Exiting as something unexpected happened."
    fi
    exit 1;
}

checkStatus(){
    if [ $? -ne 0 ]; then
	echo
	if [ $# -ge 1 ];then
	        if [ -z $2 ];then
		    errorState $1
		        else
		    errorState $1 $2
		        fi
		else
	        errorState "something unexpected happened."
		fi
    fi
}

## processing
echo "Starting the transfer"
echo "The original file is $1"
echo "The output file is $2"
echo `date`

## get the file size
origsize=`stat -c %s $1`
if [ -z $origsize ];then
    errorState "the filesize could not be determined"
fi
echo "The size of the file to transfer is: $origsize"

## check if the destination file exists
if [ -f $2 ];then
    echo "The destination file already exists."
    finalsize=`stat -c %s $2`
    echo "The destination file size is: $finalsize"
    if [ $finalsize -lt $origsize ];then
	echo "Its size is smaller than expected, removing."
	rm $2
    else
	origmd5=`md5sum $1 | cut -d' ' -f1`
	echo "The original md5 is: $origmd5"
	finalmd5=`md5sum $2 | cut -d' ' -f1`
	echo "The destination md5 is: $finalmd5"
	if [ "$finalmd5" != "$origmd5" ];then
	    echo "Its md5 is not expected, removing."
	    rm $2
	else
	    echo "Files are identical. Skipping transfer."
	    echo "Removing the original file."
	    rm $1
	    echo "Done"
	    echo `date`
	    exit 0
	fi
    fi
fi


## transfer
cp $1 $2
checkStatus "the file could not be transfered." $2
echo "The file was transfered"

## verify the transfered file size
finalsize=$(stat -c %s $2)
if [ $finalsize -ne $origsize ];then
    errorState "the file was not transfered successfully" $2
fi
echo "The transfered file has the expected size"
echo "Removing the original file."
rm $1
echo "Done"
echo `date`

