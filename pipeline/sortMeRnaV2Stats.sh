grep "%" *.log | grep passing | awk -F_ '{print $4"_"$5,$8}' | awk '{print $1,$10}' | sed 's:[(,)]::g'
grep "%" *.log | grep -v passing | grep -v failing | awk -F_ '{print $4"_"$5,$9}' | awk '{print $1,$3}' > file.txt
R dat <- read.delim(sep=" ",file="file.txt",stringsAsFactors = FALSE,header=FALSE)
do.call(rbind,split(dat[,2],dat[,1]))
