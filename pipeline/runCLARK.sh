#!/bin/bash -l
#SBATCH --mem=300G
#SBATCH -c 32
#SBATCH -w watson

set -e
set -x

module load bioinfo-tools clark

usage(){
    echo "SE mode: $0 taxonomy <infiles> <outfiles>"
    echo "PE mode: $0 taxonomy <fw_files> <rv_files> <outfiles>"
    echo
    echo "<infiles>, <outfiles>, <fw_files> and <rv_files>"
    echo "MUST be files that specify the PATH to GZipped"
    echo "files line by line. CLARK will create an output"
    echo "file for line 1 of <outfiles> using line 1 of"
    echo "<infiles> or the mates in line 1 of <fw_files>"
    echo "and <rv_files>. You MUST supply these as absolute"
    echo "paths."
    exit 1
}

# Sanity checks:

REFDIR="/mnt/picea/storage/reference/clark"

# Check argument number (must be 4)
if [ $# -ne 4 ] && [ $# -ne 3 ]
then
    usage
fi

# Check if data is paired or single end
if [ $# -eq 3 ]
then
    MODE="SE"
else
    MODE="PE"
fi

# Set up CMD arguments
if [ $MODE = "PE" ]; then
    taxonomy=$1
    fw_reads=$2
    rv_reads=$3
    out_file=$4
else
    taxonomy=$1
    fw_reads=$2
    out_file=$3
fi

# Check that all files exist (fw, [rv], out)
if [ ! -f $fw_reads ]
then
    echo "File $fw_reads does not exist."
    exit 1
fi

if [ ! -f $rv_reads ] && [ $MODE = "PE" ]
then
    echo "File $rv_reads does not exist."
    exit 1
fi

if [ ! -f $out_file ]
then
    echo "File $out_file does not exist."
    exit 1
fi

if [ ! ${fw_reads:0:1} = "/" ]; then
    fw_reads=$(pwd)/$fw_reads
fi

if [ ! ${out_file:0:1} = "/" ]; then
    out_file=$(pwd)/$out_file
fi

if [ ! ${rv_reads:0:1} = "/" ] && [ $MODE = "PE" ]; then
    rv_reads=$(pwd)/$rv_reads
fi

if [ ! $taxonomy = "genus" ] && [ ! $taxonomy = "species" ]
then
    echo "Taxonomy has to be either 'genus' or 'species'."
    exit 1
fi

# If the out directory does not exist, attempt to create it
cat $out_file | while read line
do
    out_dir=$(dirname $line)
    if [ ! ${out_dir:0:1} = "/" ]; then
	echo "Please only supply absolute paths in all files."
	exit 1
    fi
    mkdir -p $out_dir
done

# Check that all the paths are absolute and files exist
cat $fw_reads | while read line
do
    if [ ! -f $line ]; then
	echo "File: $line in $fw_files is not a valid file."
	exit 1
    fi
    if [ ! ${line:0:1} = "/" ]; then
	echo "File: $line in $fw_files is not an absolute path."
	echo "Please only supply absolute paths in all files."
	exit 1
    fi
done

if [ $MODE = "PE" ]; then
    cat $rv_reads | while read line
    do
	if [ ! -f $line ]; then
	    echo "File: $line in $rv_files is not a valid file"
	    exit 1
	fi
	if [ ! ${line:0:1} = "/" ]; then
	    echo "File: $line in $rv_files is not an absolute path"
	    echo "Please only supply absolute paths in all files."
	    exit 1
	fi
    done
fi

# #################################################
# FROM HERE ON OUT WE NEED TO CLEAN ON EXIT
# #################################################

clean()
{
    if [ -d $tempdir ]; then
	rm -rf $tempdir
    fi
}

trap clean SIGINT SIGHUP SIGTERM EXIT

# Set up a working directory
HASH=$(echo $(date +"%m%d%y%H%M%S%N")$RANDOM | md5sum | cut -d ' ' -f 1)
tempdir="/tmp/clark${HASH:0:10}"
# Make sure it really doesn't exist
while [ -f $tempdir ]
do
    HASH=$($(date +"%m%d%y%H%M%S%N")$RANDOM | md5sum | cut -d ' ' -f 1)
    tempdir="/tmp/clark${HASH:0:10}"
done
mkdir $tempdir

# Copy clark executables to current directory
bindir=$(dirname $(which classify_metagenome.sh))
cp -r $bindir/* $tempdir
cd $tempdir

# Create temporary files for Fw reads (pipes), Rv reads (pipes)
cat $fw_reads | while read line
do
    tmp=$(tempfile -d $tempdir)
    gzip -d -c $line > $tmp
    echo $tmp >> $tempdir/my_fw_files.txt
done

if [ $MODE = "PE" ]; then
    cat $rv_reads | while read line
    do
	tmp=$(tempfile -d $tempdir)
	gzip -d -c $line > $tmp
	echo $tmp >> $tempdir/my_rv_files.txt
    done
fi

# Set targets
$tempdir/set_targets.sh $REFDIR custom "--$taxonomy"

# Classify metagenome
if [ $MODE = "SE" ]; then
    $tempdir/classify_metagenome.sh --spaced -P\
				    $tempdir/my_fw_files.txt $tempdir/my_rv_files.txt\
				    -n 32 -k 31 -R $out_file
else
    $tempdir/classify_metagenome.sh --spaced -O\
				    $tempdir/my_fw_files.txt\
				    -n 32 -k 31 -R $out_file
fi
