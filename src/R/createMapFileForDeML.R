#' ---
#' title: "Reverse-complementing the index I1 and export for deML"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Setup
#' Set the working dir
setwd(".")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir=".")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(Biostrings))

#' # Input
spec <- matrix(c(
  "help",  "h" , 0, "logical"  , "display a help message",
  "chunks",  "c" , 255, "integer"  , "the number of sample per map",
  "outdir", "o", 1, "character","the out dir",
  "file", "f", 1, "character","a map file",
  "revComp", "r", 0, "logical", "reverse complement Index1"
), byrow=TRUE, ncol=5)

# get the options
opt <- getopt(spec)

# local functions
".abort" <- function(msg){
  message(paste("ERROR:",msg))
  .help(1)
}

".help" <- function(status=0){
  message(getopt(spec, usage=TRUE))
  message("\nThe -f option expects a map file with three columns")
  message("\ni.e. #Index1, Index2 and Name, tab separated.")
  quit(save="no",status=status);
}

# help
if ( !is.null(opt$help) ) {
  .help()
}

# validation
if(is.null(opt$file) | is.null(opt$outdir)){
  .abort("The -f and -d arguments are required.")
}

if(!file.exists(opt$file)){
  .abort("The -f options should be an existing file")
}

if(!file.exists(opt$outdir)){
  .abort("The -d options should be an existing directory")
}

if ( is.null(opt$chunks) ) {
  opt$chunks=255
}

message(opt$outdir)
message(opt$file)
message(opt$chunks)

dat <- read.delim(opt$file,header=TRUE,as.is=TRUE)

if(!is.null(opt$revComp)){
  dat[,1] <- as.character(reverseComplement(DNAStringSet(dat[,1])))
}

if(any(duplicated(dat[,1:2]))){
  .abort("Duplicated index pairs identified")
}

# write in chunks of 255 by default
# deML opens 8 filehandles per sample (plus a few extra) and system limit is commonly 2048
print(opt$chunks)
chks <- breakInChunks(nrow(dat),opt$chunks)

dev.null <- sapply(chks,function(p,d,o,f){
  f <- file.path(o,paste0(sub("\\..*$","",f),"_",min(p),"-",max(p),"_deML.txt"))
  write("#Index1\tIndex2\tName",file=f)
  write.table(d[p,],file=f,append = TRUE,sep="\t",
              row.names = FALSE, quote = FALSE,
              col.names = FALSE)
},dat,opt$outdir,basename(opt$file))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
