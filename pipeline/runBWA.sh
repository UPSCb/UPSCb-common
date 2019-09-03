#! /bin/bash -l
#SBATCH -A b2010042
#SBATCH -p node -n 16
## 2 days for fosmidpools, 4 for diploid
#SBATCH -t 2-00:00:00
#SBATCH -J diginormBWAPE
## usw fat for fosmidpools and masterassembly, 72 for the diploid
##SBATCH -C fat
##SBATCH -C mem72GB
#SBATCH --mail-user nicolas.delhomme@plantphys.umu.se
#SBATCH --mail-type=ALL

#if [ ! -z $SLURM_SUBMIT_DIR ]; then
    module load bioinfo-tools bwa samtools
#fi

GZIP=
PROC=16
N=0
O=0

## a sub 
usage(){
echo >&2 \
"
	  usage: $0 [options] <fwd fq> <rv fq> <genome> <out>
	  options:
	     -g: set the path to a tmp directory if the read files are compressed
             -n: max #diff (int) or missing prob under 0.02 err rate (float) [$N]
             -o: maximum number or fraction of gap opens [$O]
             -t: the number of thread to use. Default to $PROC
          note: All other options are set to default. Ask Nicolas Delhomme to extend the present script if you have other requirements.
"
exit 1
}

## get the mode
while getopts g:n:o:t: opt
do
    case "$opt" in
	g) GZIP=$OPTARG;;
	n) N="$OPTARG";;
	o) O="$OPTARG";;
	t) PROC="$OPTARG";;
	\?)		# unknown flag
      	    usage;;
    esac
done
shift `expr $OPTIND - 1`

if [ $# != 4 ]; then
    echo "This function takes 2 fastq files, a genome directory prefix, and an ouput directory"
    usage
fi

if [ ! -f $1 ]; then
    echo "The first argument should be the fwd fastq file"
    usage
fi
left=$1

if [ ! -f $2 ]; then
    echo "The second argument should be the reverse fastq file"
    usage
fi
right=$2

if [ ! -f $3.bwt ]; then
    echo "The third argument should be the BWA genome index prefix"
    usage
fi
genome=$3

if [ ! -d $4 ]; then
    echo "The forth argument should be the output directory"
    usage
fi
out=$4

if [ ! -z $GZIP ]; then
    if [ ! -d $GZIP ]; then
	echo "The -g option expects a tmp directory as argument"
	usage
    fi
    gunzip -c $left > $GZIP/`basename ${left//.gz/}`
    left=$GZIP/`basename ${left//.gz/}`
    gunzip -c $right > $GZIP/`basename ${right//.gz/}`
    right=$GZIP/`basename ${right//.gz/}`
fi


fnam=`basename ${left//.f*q/}`

## now run
bwa aln -o $O -n $N -t $PROC $genome $left > $out/${fnam}-1.sai 2> $out/${fnam}-1.err

bwa aln -o $O -n $N -t $PROC $genome $right > $out/${fnam}-2.sai 2> $out/${fnam}-2.err

bwa sampe $genome $out/${fnam}-1.sai $out/${fnam}-2.sai $left $right | samtools view -bS - | samtools sort - $out/${fnam} 2> $out/${fnam}-bwa.err

samtools index $out/$fnam.bam

## clean up
if [ ! -z $GZIP ]; then
rm $left $right
fi

## rm all.DigiNorm.C20_k20_N4_x2e9-1.sai
## rm all.DigiNorm.C20_k20_N4_x2e9-2.sai
