## load the libraries
library(Rsamtools)
library(GenomicAlignments)

setwd("/mnt/picea/projects/aspseq/sex/STAR/Potra01")

bamFileList <- BamFileList(dir(".",pattern="*_STAR.bam$")[1:2],yieldSize = 10^6)

system.time( gList <- mclapply(bamFileList,function(bamFile){
  
  ## open the bamFile
  open(bamFile)
  
  ## create an empty container
  gPos <- GRanges()
  
  ## fill it by processing 1M reads at a time
  while(sum(elementLengths(records <- scanBam(bamFile,
                                              param=ScanBamParam(
                                                flag=scanBamFlag(isUnmappedQuery=FALSE),
                                                what=scanBamWhat()[c(3,5,8)]
                                              ))[[1]]
  ))>0){
    message("Processing 1M reads")
    message(paste("Reading from record",records$rname[1],records$pos[1],records$cigar[1]))
    nOps <- cigarRangesAlongReferenceSpace(records$cigar,ops="N")
    sel <- elementLengths(nOps) != 0
    gPos <- c(gPos,GRanges(seqnames=rep(records$rname[sel],elementLengths(nOps)[sel]),
                           ranges=unlist(shift(nOps[sel],records$pos[sel]))))
  }
  
  ## close the bamFile
  close(bamFile)
  
  ## return all gap positions
  return(gPos)
},mc.cores=2L))

## get all unique positions across all samples
uPos <- unique(Reduce("c",gList))

## and it's a small object
print(object.size(uPos),unit="Mb")

## Finally create the "count table"
## Because the genome used as >200,000 scaffolds, the summarizeOverlaps 
## function does not perform well, so we revert to a quick and dirty fix
ref <- paste(seqnames(uPos),start(uPos),end(uPos),sep="-")

system.time(count.table <- do.call("cbind",mclapply(gList,function(grng,ref){
  tab <- table(paste(seqnames(grng),start(grng),end(grng),sep="-"))
  tab[match(ref,names(tab))]
},ref,mc.cores=2L)))

rownames(count.table) <- ref
head(count.table)

## sessionInfo()
sessionInfo()
