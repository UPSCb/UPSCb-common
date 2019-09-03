#!/bin/bash

## a sub 
usage(){
echo >&2 \
		    "
For --semap or --pereport:
	  usage: $0 --semap [options] <fwd fq> <rv fq> <genome dir> <genome> <gene model>
	  usage: $0 --pereport [options] <bam file> <genome dir> <genome> <gene model>
	  options:
	     -d: the output directory. Default to the current directory
	     -f: a file name for the ini file. By default a temporary file is created and erased at the end of the process	     
	     -g: set if the read files are compressed
             -m: mock - do not run
	     -n: the output file name prefix. Default to 'semap[PID]'
	     -p: whether to use the Paired-End mode or not. Defaut to false.
             -t: the number of thread to use. Default to 8
          note: All other options are set to default. Ask Nicolas Delhomme to extend the present script if you have other requirements.

For anything else, check the standard FusionMap command
"
exit 1
}

gzip="False"
paired="False"
ini=
mock="False"
odir=`pwd`
oname="semap$$"
thread=8

## get the mode
mode=$1
format="FASTQ"
if [ "$mode" == "--semap" -o "$mode" == "--pereport" ]; then
    shift
    
    while getopts gpf:d:mn:t: opt
    do
	case "$opt" in
	    g) gzip="True";;
	    p) paired="True";;
	    f) ini="$OPTARG";;
	    d) odir="$OPTARG";;
	    m) mock="True";;
	    n) oname="$OPTARG";;	    
	    t) thread="$OPTARG";;
	    \?)		# unknown flag
      		usage;;
	esac
    done
    shift `expr $OPTIND - 1`

    ## sanity
    case "$mode" in 
	--semap)
	    if [ ! -f $1 ]; then
		echo "The forward fastq file cannot be found."
		usage
	    fi
	    if [ ! -f $2 ]; then
		echo "The reverse fastq file cannot be found."
		usage
	    fi
	    if [ ! -d $3 ]; then
		echo "Error: the third argument is not a directory (--semap mode)."
		usage
	    fi;;	
	--pereport)
	    if [ ! -f $1 ]; then
		echo "The sam file cannot be found."
		usage
	    fi
	    if [ ! -d $2 ]; then
		echo "Error: the second argument is not a directory (--pereport mode)."
		usage
	    fi
	    format="BAM"
	    ;;
    esac

## create an ini file
    if [ -z "$ini" ]; then
	ini=`mktemp`
    fi
    
    echo "<Files>" > $ini
    echo $1 >> $ini
    shift
    if [ "$mode" == "--semap" ]; then
	echo $1 >> $ini
	shift
    fi
    echo "" >> $ini
    echo "<Options>" >> $ini
## Possible values: True, False. Default value=True
    echo "RnaMode=True" >> $ini
## Possible values: True, False. Default value=False
    echo "PairedEnd=$paired" >> $ini
## Possible values: depends on your machine. Default value=8 
    echo "ThreadNumber=$thread" >> $ini
## Possible: True, False. Default=False, read through events<50Kb may missed if this option is on.
    echo "SearchNovelExonJunction=True" >> $ini
## Possible values: FASTQ, QSEQ, FASTA. Default value=FASTQ for semap; SAM for pereport (we use BAM)
    echo "FileFormat=$format" >> $ini
## Possible values: True, False. Default value=True
    echo "AutoPenalty=True" >> $ini
## Possible values: 0-100. Default value=2
    echo "FixedPenalty=2" >> $ini 
## Possible values: 0-100. Default value=2
    echo "IndelPenalty=2" >> $ini 
## Possible values: True, False. Default value=False
    echo "DetectIndels=False " >> $ini 
## Possible values: True, False. Default value = True for 64-bit OS and False for 32-bit OS
    echo "Use32BitMode=False" >> $ini 
## Possible values: 1-50. Default value=10
    echo "MaxMiddleInsertionSize=10" >> $ini 
## Possible values: 1-500. Default value=10
    echo "MaxMiddleDeletionSize=10" >> $ini 
## Possible values: 1-50. Default value=10
    echo "MaxEndInsertionSize=10" >> $ini 
## Possible values: 1-500. Default value=10
    echo "MaxEndDeletionSize=10" >> $ini 
## Possible values: 3-10. Default value=3
    echo "MaxDistalEndSize=3" >> $ini 
## Possible values: True, False. Default value=True
    echo "TrimByQuality=False" >> $ini 
## Possible values: 20-65536. Default value = 1024
    echo "ReadTrimSize=1024" >> $ini 
## Possible values: 0-20. Default value=2
    echo "ReadTrimQuality=2" >> $ini 
## Possible values: Automatic, Illumina, Sanger. Default value=Automatic
    echo "QualityEncoding=Automatic" >> $ini 
## Possible values: True, False. Default value=False
    echo "Gzip=$gzip" >> $ini 
## Possible values: 15-50. Default value=25
    echo "MinimalFusionAlignmentLength=25" >> $ini 
## Possible values: True, False. Default Value=True
    echo "FilterUnlikelyFusionReads=False" >> $ini 
## Possible values:0-500000. Default value=5000
    echo "MinimalFusionSpan=50000" >> $ini 
## Possible values: 1-10000, Default value =2
    echo "MinimalHit=2" >> $ini 
## Minimal rescued read number. Default value = 1
    echo "MinimalRescuedReadNumber=1" >> $ini 
## Possible values: True, False. Default value = False
    echo "ReportUnannotatedFusion=True" >> $ini 
## Possible values: True, False. Default value = True
    echo "OutputFusionReads=True" >> $ini 
## Possible values: 1-1000. Default value=1
    echo "FusionReportCutoff=1" >> $ini 
## Possible values: 01-. Default value = 2
    echo "NonCanonicalSpliceJunctionPenalty=2" >> $ini
## Possible values: “DefaultList”, “CustomizedList”, or “None”
    echo "FilterBy=None" >> $ini
##
## DefaultFilterListVersion=v1 
## Specify the file path if FilterBy=CustomizedList
## FilterGeneListFileName= 
## Specify the file path if FilterBy=CustomizedList
## FilterGeneFamilyFileName=

    echo "" >> $ini
    echo "<Output>" >> $ini
    echo "OutputName=$oname" >> $ini
    echo "OutputPath=$odir" >> $ini
    
## start the script
    if [ "$mock" == "False" ]; then
	mono /home/delhomme/opt/FusionMap_2014-01-01/bin/FusionMap.exe --semap $@ $ini
    else
	echo "Mock run finished. The ini file is: $ini"
    fi
else
## otherwise go on
## start the script
    mono /home/delhomme/opt/FusionMap_2014-01-01/bin/FusionMap.exe $@
fi
