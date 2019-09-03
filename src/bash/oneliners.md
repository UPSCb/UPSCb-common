# Useful bash one-liners

Edit this file to add oneliners (you can cheat with `;` ) you find useful. If there is no appropriate sub-heading, feel free to create one.

# Table of contents

1. [SLURM related](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#slurm-related)
    * [Cancel all SLURM jobs that you have currently running](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#cancel-all-slurm-jobs-that-you-have-currently-running)
    * [Delete all SLURM logs in the current directory (and all sub-directories)](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#delete-all-slurm-logs-in-the-current-directory-and-all-sub-directories)
2. [File and folder management](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#file-and-folder-management)
    * [Find all human-readable files with pattern in text](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#find-all-human-readble-files-with-pattern-in-text)
3. [Miscellanea](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#miscellanea)
4. [Docker](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#docker)
5. [RNA-Seq Preprocessing](https://microasp.upsc.se/root/UPSCb/blob/master/src/bash/oneliners.md#rna-seq-preprocessing)

## SLURM related

### Cancel all SLURM jobs that you have currently running

```bash
squeue -h -u $USER | awk '{ print $1 }' | xargs -I {} scancel {}
```

### Delete all SLURM logs in the current directory (and all sub-directories)

```bash
find . -type f -name "slurm-*.out" -user $USER -exec rm {} \;
```

### Re-assign the requested node

```bash
squeue -h -u delhomme | awk '{if($5=="PD"){print $1}}' | xargs -I {} sudo scontrol update JobId={} ReqNodeList=picea
```

## File and folder management

### Find all human-readble files with pattern in text

```bash
# Note the pattern that should be searched at the end of the command.
find . -type f -exec sh -c 'file -b {} | grep text &>/dev/null' \; -print | xargs -I {} grep -H "PATTERN" {}
```

## Miscellanea

### Biology related

#### Counting all bp in a fasta file

```bash
grep -v ">" CCS_FrameDP.fa | sed ':a;N;$!ba;s/\n//g' | wc
```

#### Finding the intron size range

```bash
grep intron TAIR10_GFF3_genes_transposons_introns.gff | awk 'BEGIN{max=0;min=0}{siz=$5-$4; if(siz>max){max=siz}; if(siz<min){min=siz}} END {print min+1 "-" max+1}'
```

#### Getting the number of unmapped reads from STAR

In a _star_ directory or by changing the _find_ command to point to a _star_
directory, it counts the number of unmapped reads (no matter why) per sample.

```bash
for f in `find . -name *final.out`; do echo `basename ${f//Log.final.out/}`; grep umber $f | grep reads | cut -f2 | awk '{if(NR==1){tot=$1}else{tot-=$1}}END{print tot}'; done
```

#### convert sam to bam
```bash
find . -name "*.sam" | xargs -I {} srun bash -c 'samtools view -b $0 | samtools sort - $0' {}
```

### Computer Science related

#### Finding files longer than 80bp

```bash
awk 'length > 80 {print FILENAME "(" FNR "): " $0}'
```

#### Replacing TAB with 4 spaces

```bash
find DESCRIPTION NAMESPACE NEWS vignettes/RnaSeqTutorial.Rnw man/*.Rd R/*.R | xargs -I {} bash -c 'perl -p -e "s/\t/    /g" $0 > $0.mod' {}
```

## Docker

### Clean all stopped containers:

```bash
docker rm $(docker ps -a -f "status=exited" -q)
```

### Clean all dangling images:

```bash
docker rmi $(docker images -f "dangling=true" -q)
```

## RNA-Seq preprocessing

### Count the number of reads after the pipeline run
```bash
cd fastqc
find raw/multiview/*_1_fastqc -name fastqc_data.txt | sort | xargs -I {} bash -c 'echo $0 | awk -F_ "{printf \"%s_%s \",\$4,\$5}" ; grep "Total Sequences" $0 | awk "{print \$3}"' {}
```