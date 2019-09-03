find $1 -name "fastqc_data.txt" | xargs -I {} bash -c 'echo $0 $(grep "Total Sequences" $0)' {}

#find . -name "fastqc_data.txt" | xargs -I {} bash -c 'echo $0 $(grep "Total Sequences" $0)' {} | awk '{id=$1;gsub("\\./|_[1,2]_fastqc.*","",id);print id, $4}' | sort | uniq