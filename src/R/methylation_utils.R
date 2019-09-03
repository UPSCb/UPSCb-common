combineContexts <- function(in.prefix){
  # Combine several methylation call files from methylkit
  # into a sorted BED-like file
  require(data.table)
  CpGfile <- paste(in.prefix,"_CpG.txt",sep="")
  CHGfile <- paste(in.prefix,"_CHG.txt",sep="")
  CHHfile <- paste(in.prefix,"_CHH.txt",sep="")
  if(any( ! sapply(c(CpGfile,CHHfile,CHGfile), file.exists))) {
    stop(paste("File error. Please make sure your prefix is corrent and your",
               "file endings match '_CpG|CHH|CHG\\.txt'"),sep="\n")
  }
  CpG <- fread(CpGfile,showProgress = FALSE)
  CHG <- fread(CHGfile,showProgress = FALSE)
  CHH <- fread(CHHfile,showProgress = FALSE)
  CpG$context <- "CpG"
  CHG$context <- "CHG"
  CHH$context <- "CHH"
  comb <- rbind(CpG,CHG,CHH)
  res <- comb[, .(chr,start=base,end=base,coverage,
                  intT = round(coverage*freqT/100),context)]
  setkey(res,chr,start,end,context)
  return(res)
}

mergeMKit <- function(...){
  # Merge objects returned by combineContexts as a left outer join and retain only
  # site present in all samples
  inputs <- list(...)
  if( length(inputs) < 2 ){
    inputs <- as.list(unlist(inputs))
    inputs.names <- unlist(inputs)
    if( length(inputs) < 2 ){
      stop("Not sure how to merge this")
    }
  } else {
    inputs.names <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
  }

  l <- length(inputs)
  res <- data.table()
  if(class(inputs[[1]])=="character"){
    message("Reading files")
    inputs <- lapply(seq_along(inputs),function(i){
      inp <- inputs[[i]]
      message(paste("Read",i,"file(s)"))
      return(combineContexts(inp))
    })
  }
  for(f in 1:l){
    message(paste("Merge step",f,"of",l))
    if(f==1){
      res <- inputs[[f]]
    } else {
      res <- merge(res,inputs[[f]])
      sample.cols <- unlist(lapply(inputs.names[1:f],
                                   paste,c(".coverage",".T"),sep=""))
      setnames(res,c("chr","start","end","context",sample.cols))
    }
  }
  sample.cols <- unlist(lapply(inputs.names[1:f],
                               paste,c(".coverage",".T"),sep=""))
  setnames(res,c("chr","start","end","context",sample.cols))
  setkey(res,chr,start,end,context)
  return(res)
}