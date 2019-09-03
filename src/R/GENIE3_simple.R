tic <- Sys.time()

difft <- function(){
  toc <- Sys.time()
  ti <- round(difftime(toc,tic, units="secs"))
  return(ti)
}

statusUpd <- function(mess){
  message(paste("[\t",difft(),"]\t",mess))
}

statusUpd("Loading packages")

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(igraph)
  library(parallel)
  library(matrixStats)
})

statusUpd("Getting arguments")

args <- commandArgs(trailingOnly = TRUE)

inf <- args[1]
outf <- args[2]
cores <- args[3]

statusUpd("Checking file permissions")
if(! file.create(outf)){
  stop(paste("Error creating output file.",
             "Check that you have the appropriate permissions"))
} else {
  file.remove(outf)
}

statusUpd("Testing multicore affinity")
if(length(mcaffinity()) < cores){
  stop(paste("Length of mc affinity < number of reserved cores.",
             "Did you reserve enough SLURM resources?"))
} else {
  mcaffinity(seq(detectCores()))
}

statusUpd("Reading input data")
dat <- fread(inf)
samples <- nrow(dat)
dat <- dat[,colSums(sign(as.matrix(dat))) > sqrt(samples),with=FALSE]
genes <- colnames(dat)

outd <- dirname(outf)
tempd <- paste(outd,"/","tmp",sep="")
dir.create(tempd)

statusUpd("Starting rf computation")
res <- mclapply(seq(1,ncol(dat)), mc.cores = cores,
                function(gene){
                  if(gene%%100==0){
                    perc <- round(gene/length(genes)*100,2)
                    statusUpd(paste("Completed ",gene," genes. ",
                                    perc,"%",sep=""))
                  }
                  t <- tempfile(tmpdir = tempd)
                  i <- gene
                  s <- setdiff(seq(1,ncol(dat)),i)
                  y <- unlist(dat[,i,with=FALSE])
                  x <- dat[,s,with=FALSE]
                  im <- importance(randomForest(x,y,importance=TRUE))[,2]/samples
                  dt <- data.table("i" = genes[i],"j" = genes[s],
                                   "weight" = im)
                  dt <- dt[weight>0]
                  write.table(dt,t,quote = FALSE,col.names=FALSE,
                              row.names = FALSE, sep="\t")
                  return(TRUE)
                })
