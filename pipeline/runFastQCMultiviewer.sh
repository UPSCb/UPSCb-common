#!/bin/bash

set -ex

usage(){
    echo >&2 \
"
    runFastQCMultiviewer.sh <fastqc dir>

    Arguments:
        fastqc dir - the directory containing the FastQC reports

    Note:
        The UPSCb Environment Variable needs to be set to your
        Git UPSCb checkout dir.

    Details:
        It unzip the fastqc files into a directory called multiview
        and create a multiview.html file.
"
exit 1
}

if [ $# -ne 1 ]; then
    echo "This function takes one argument"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be the directory containing the FastQC reports"
    usage
fi

py=${SLURM_SUBMIT_DIR:-$(dirname $0)}/../src/python/fastQCmultiviewer.py

if [ ! -f $py ]; then
    echo "Fixme; the .. part of the path, also below..."
    usage
fi

mkdir -p $1/multiview

find $1 -name "*.zip" -type f -exec unzip -f -d $1/multiview "{}" \;

python $py -out_file $1/multiview.html -in_dir $1/multiview

