#!/bin/bash
/proj/uppstore2017112/b2017108_nobackup

# count file extensions
find GenomicSelection -type f ! -name "*.gz" ! -name "*.pdf" ! -name "*.idx" | xargs -I {} bash -c 'echo ${0##*.}' {} | sort | uniq -c

# cleanup
find . -type l -delete
find . -name "._*" -delete
find . -type f -name "slurm*.out" -delete
find . -name ".DS_Store" -delete
find . -name ".RData*" -delete
find . -name ".nfs*" -delete
find . -name "*.bak" -delete
find . -name "*.save*" -delete
find . -name "*.swp*" -delete
find . -name "*core.[0-9]*" -delete


# compression extension based
find /proj/uppstore2017112/b2017108_nobackup/GenomicSelection -type f -name "*.vcf" > vcf.list
find /proj/uppstore2017112/b2017108_nobackup/GenomicSelection -type f -name "*.txt" > text.list
find /proj/uppstore2017112/b2017108_nobackup/GenomicSelection -type f -name "*.out" >> text.list 
find /proj/uppstore2017112/b2017108_nobackup/GenomicSelection -type f -name "*.mendel" > mendel.list
find /proj/uppstore2017112/b2017108_nobackup/GenomicSelection -type f -name "*.fa" -o -name "*.fam" > fasta.list
find $(realpath .) -name "*bed" -o -name "*.csv" -o -name "*.dat" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.gff" -o -name "*.gff3" -o -name "*.gtf" -o -name "*.indv" -o -name "*.intervals" -o -name "*.junction" -o -name "*.log" -o -name "*.mate1" -o -name "*.mate2" -o -name "*.nosex" -o -name "*.pos" -o -name "*.sam" -o -name "*.sample" -o -name "*.snp" -o -name "*.tab" -o -name "*.text" -o -name "*.tranches" -o -name "*.tsv" -o -name "*.xml" > mix.list

# split large files
split -l 1000 vcf.list vcf
split -l 1000 text.list txt
split -l 1000 mix.list mix

# run as job arrays
sbatch -A snic2019-8-261 -t 6:00:00 -a 0-48 -o fa.out -e fa.err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh fasta.list
sbatch -A snic2019-8-261 -t 6:00:00 -a 0-32 -o mendel.out -e mendel.err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh mendel.list
find . -maxdepth 1 -name "vcfa[a-f]" -! -name "vcfa.*" -exec sbatch -A snic2019-8-261 -t 6:00:00 -a 0-999 -o "{}".out -e "{}".err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh "{}" \;
find . -maxdepth 1 -name "vcfa[g-j]" -! -name "vcfa.*" -exec sbatch -A snic2019-8-163 -t 6:00:00 -a 0-999 -o "{}".out -e "{}".err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh "{}" \;
find . -maxdepth 1 -name "txta[a-c]" -! -name "txta.*" -exec sbatch -A snic2019-8-163 -t 6:00:00 -a 0-999 -o "{}".out -e "{}".err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh "{}" \;
find . -maxdepth 1 -name "txta[d-i]" -! -name "txta.*" -exec sbatch -A snic2019-8-124 -t 6:00:00 -a 0-999 -o "{}".out -e "{}".err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh "{}" \;
find . -maxdepth 1 -name "mixa[a-h]" -! -name "mixa.*" -exec sbatch -A snic2019-8-124 -t 6:00:00 -a 0-999 -o "{}".out -e "{}".err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh "{}" \;

# xxhash
find . -size +1G | xargs -P 10 -I {} ~/bin/xxhsum {} > xxhash.list &

# manual cleanup of duplicated files

# compression file-type based
# this command finds file that are non binary (- capital i, lower case L); the {} + is a faster paradgim to "{}" \; : find . -type f -exec grep -Il . {} + 
find . -type f ! -name "*.gz" -size +100M -exec grep -Il . {} + > ascii.list
sbatch -A snic2019-8-124 -t 6:00:00 -a 0-493 -o ascii.out -e ascii.err ~/Git/UPSCb/pipeline/runAsArray.sh ~/Git/UPSCb/pipeline/runGzip.sh ascii.list 






