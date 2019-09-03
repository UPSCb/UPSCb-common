library(stringr)
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  cat("Usage: Rscript updateVcf.R <conversionTable> <infile> <outfile>\n")
  stop("Usage error")
}

conv <- args[1]
infile <- args[2]
outfile <- args[3]
cat("Reading vcf....\n")
f <- readLines(infile)
cat("Reading conversion table....\n")
conv <- read.table(conv, header = FALSE, stringsAsFactors = FALSE)

cat("Changing vcf header....\n")

reheader <- function(f){
  contigIndex <- grepl("^\\#\\#contig=<ID=.*,",f)
  header <- gsub(",","",gsub("##contig=<ID=","",str_extract(f[contigIndex],"^\\#\\#contig=<ID=.*,")))
  header <- conv[,2][match(header,conv[,1])]
  newHeader <- paste('##contig=<ID=',header,',',sep = '')
  repl <- str_replace(f[contigIndex],"^\\#\\#contig=<ID=.*,",newHeader)
  f[contigIndex] <- repl
  return(f)
}

f <- reheader(f)

cat("Changing rest of vcf....\n")

retail <- function(f){
  tailindex <- !grepl("\\#",f)
  newtail <- conv[,2][match(gsub("\\t","",str_extract(f[tailindex],".*?\\t")),conv[,1])]
  newtail <- paste(newtail,"\t",sep="")
  repl <- str_replace(f[tailindex],".*?\\t",newtail)
  f[tailindex] <- repl
  return(f)
}

f <- retail(f)

cat("writing to disk",outfile,"....\n")

write(x = f,file = outfile)
