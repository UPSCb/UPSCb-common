## load the optparse library
suppressPackageStartupMessages(library(optparse))

### ================ functions
getIDs <- function(fl){
  library(ShortRead)
  fs <- FastqStreamer(fl,n=1e6)
  ids <- BStringSet()
  while(length(fq <- yield(fs))){
      ids <-c(ids,BStringSet(sub("/[1,2].*","",id(fq))))
      ## strange ids sometimes... 
    ## ids <- c(ids,narrow(id(fq),start=1,end=width(id(fq))-2))
  }
  return(ids)
}

## because it's slow, only write 1,000,000 reads at once
writeFqs <- function(i,ins,cid,pouts,souts){
  library(ShortRead)
  fs <- FastqStreamer(ins[[i]],n=1e6)
  while(length(fq <- yield(fs))){
      sel <- is.na(match(BStringSet(sub("\\#.*","",id(fq))),cid))
      ## paired
      writeFastq(fq[!sel],file=pouts[[i]],mode="a")
      ## single
      writeFastq(fq[sel],file=souts[[i]],mode="a")
      ##    writeFastq(fq[!is.na(match(narrow(id(fq),start=1,end=width(id(fq))-2),cid))],file=outs[[i]],mode="a")  
  }
  return(TRUE)
}

### ================ main
## define the arguments
option_list <- list(
                    make_option(c("-f", "--forward"), type="character", default="",
                                help="The forward fastq file"),                    
                    make_option(c("-r", "--reverse"), type="character", default="",
                                help="The reverse fastq file"),
                    make_option(c("-n","--nnodes"), type="integer", default=2,
                                help="Define the number of processors to use."),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print extra output [default]"))

## get the options
opt <- parse_args(OptionParser(option_list=option_list))

## load the rest of the libs
library("ShortRead")
library("parallel")

if(opt$verbose){
  message("Starting")
}

## check them
stopifnot(file.exists(opt$forward))
stopifnot(file.exists(opt$reverse))
ins <- list(opt$forward,opt$reverse)

## we just expect fq and fastq extension here
pouts <- sub("1\\.","_paired_1\\.",ins)
pouts <- sub("2\\.","_paired_2\\.",pouts)
souts <- sub("paired","single",pouts)
pouts <- sub("\\.gz$","",pouts)
souts <- sub("\\.gz$","",souts)

## check that the file are not there yet
if(opt$verbose){
  message("Cleaning")
}
sapply(c(pouts,souts),function(fl){if(file.exists(fl)){file.remove(fl)}})

### three main steps

## 1. process the files by chuncks to get the IDs
if(opt$verbose){
  message("Getting the ID list")
}
cluster <- makePSOCKcluster(opt$nnodes)
ids <- clusterApply(cluster,ins,getIDs)

## 2. determine the IDs to keep
if(opt$verbose){
  message("Finding the common IDs")
}
cids <- intersect(ids[[1]],ids[[2]])
## uids <- setdiff(union(ids[[1]],ids[[2]]),cids)

## 3. write the pairs
if(opt$verbose){
  message("Writing the paired and unpaired files")
}
res <- clusterApply(cluster,list(1,2),writeFqs,ins,cids,pouts,souts)

## 4. write the unique
#if(opt$verbose){
#  message("Writing the unpaired files")
#}
#res <- clusterApply(cluster,list(1,2),writeFqs,ins,uids,souts)

## 4. stop the cluster
stopCluster(cluster)

if(opt$verbose){
  message("Done")
}

quit(save="no")

