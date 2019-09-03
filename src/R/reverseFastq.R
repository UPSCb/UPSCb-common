## libs
library("ShortRead")

## a warning
message("That script is hardcoded as h..l, so if you need it, edit or enhance :-)")

## setwd
setwd("/mnt/picea/storage/projects/07_Sd_ludwigii_Project/fastq/454/8k")

## read the fq file
setMethod(f="reverse",
          signature="ShortReadQ",
          definition=function(x,...){
  x@sread <- reverseComplement(sread(x))
  x@quality@quality <- reverse(quality(quality(x)))
  return(x)
})

writeFastq(reverse(readFastq("trim_13_1.fastq.gz")),file="rev_trim_13_1.fastq")
writeFastq(reverse(readFastq("trim_15_1.fastq.gz")),file="rev_trim_15_1.fastq")
