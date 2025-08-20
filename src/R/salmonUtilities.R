#' ---
#' title: "Salmon utilities"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions that help extract
#' various information from a [salmon](https://salmon.readthedocs.io/en/latest/) run.
#' 
#' To use it in your script, assuming you have a checkout of the repos in the current dir, do:
#' ```{r, eval=FALSE}
#' library(here)
#' source("here(UPSCb-common/src/R/salmonUtilities.R"))
#' ```
#' 
#' In the following, you will first find the S4 generics of the different 
#' functions and then their implementation, as well as some description of
#' their functionality, arguments, etc.
#' 
#' # Libraries
#' 
#' The following libraries are necessary for the correct execution of the utility.
#' As this file is meant to be _sourced_, we _require_ the presence of the 
#' libraries.
stopifnot(suppressPackageStartupMessages({
  require(dplyr)
  require(ggplot2)
  require(magrittr)
  require(parallel)
  require(readr)
  require(S4Vectors)
  require(tibble)
}))

#' # Generics
setGeneric(name="extractBinary",
           def=function(salmonDir=character(0L),
                        type=c("gc","seq","pos"),
                        ...
           ){
             standardGeneric("extractBinary")
           })

setGeneric(name="extractEquivalenceClass",
           def=function(salmonDir=character(0L),
                        mc.cores=1L,
                        ...
           ){
             standardGeneric("extractEquivalenceClass")
           })

setGeneric(name="extractUniqueCountProportion",
           def=function(salmonDir=character(0L),
                        mc.cores=1L,
                        ...
           ){
             standardGeneric("extractUniqueCountProportion")
           })

setGeneric(name="plotBias",
           def=function(bias=list(),
                        type=c("gc","seq","pos"),
                        ...
           ){
             standardGeneric("plotBias")
           })

setGeneric(name="plotEqClassQC",
           def=function(eqc=list(),
                        geneIDs=character(0L),
                        ...
           ){
             standardGeneric("plotEqClassQC")
           })

setGeneric(name="plotProportion",
           def=function(tb=tibble(),
                        tx=character(0L),
                        ...
           ){
             standardGeneric("plotProportion")
           })

#' # Implementation
#' 
#' ## Arguments
#'
#' ## Functions
#' ### extractBinary
#'  
#' This function is a wrapper that based on the desired bias type (gc, seq or pos) delegates
#' to the appropriate "worker".
setMethod(f = "extractBinary", 
          signature = c("character","character"),
          definition = function(
    salmonDir=character(0L),
    type=character(0L),
    ...
          ){
            
            type <- match.arg(type,eval(formals("extractBinary")$type))
            
            # dispatch
            .extractBinary(salmonDir,type,...)
          })

setMethod(f = "extractEquivalenceClass", 
          signature = c("character"),
          definition = function(
    salmonDir=character(0L),
    mc.cores=1L,
    ...
          ){
            .validate(dir=salmonDir)
            
            files <- list.files(salmonDir,
                                recursive=TRUE,
                                full.names=TRUE,
                                pattern="eq_classes.txt.gz")
            
            lst <- mclapply(files,function(f){
              ntx <- scan(f,what=integer(),nmax=1)
              nec <- scan(f,what=integer(),skip=1,nmax=1)
              tx <- scan(f,what=character(),skip=2,nmax=ntx,sep="\n")
              eqc <- strsplit(scan(f,what=character(),skip=ntx+2,nmax=nec,sep="\n"),"\t")
              mlist <- mapply("[",eqc,sapply(S4Vectors::elementNROWS(eqc)-1,":",2))
              list(size=as.integer(sapply(eqc,"[[",1)),
                   members=data.frame(EQC=rep(1:length(mlist),elementNROWS(mlist)),
                                 TX=tx[as.integer(unlist(mlist))+1]),
                   counts=as.integer(mapply("[[",eqc,elementNROWS(eqc))))
            },mc.cores=mc.cores)
            
            names(lst) <- basename(dirname(dirname(files)))
            return(lst)
          })

setMethod(f = "extractUniqueCountProportion", 
          signature = c("character"),
          definition = function(
    salmonDir=character(0L),
    mc.cores=1L,
    ...
          ){
            .validate(dir=salmonDir)
            
            files <- list.files(salmonDir,
                                recursive=TRUE,
                                full.names=TRUE,
                                pattern="ambig_info.tsv")
            
            tb <- do.call(bind_cols,mclapply(files,function(f){
              nam <- basename(dirname(dirname(f)))
              tb <- read_tsv(f,show_col_types=FALSE) %>% 
                mutate(X=UniqueCount/(UniqueCount+AmbigCount)) %>% 
                select(X)
              colnames(tb) <- nam
              return(tb)
            },mc.cores=mc.cores))
            
            tx <- read_tsv(file.path(dirname(dirname(files[[1]])),"quant.sf"),
                           show_col_types=FALSE) %>% 
              select("Name")
            tb %<>% mutate(Transcript=unlist(tx,use.names=FALSE)) %>% 
              relocate(Transcript,1)  
            
            return(tb)
          })

setMethod(f = "plotBias",
          signature = c("list","character"),
          definition = function(bias=list(),
                                type=character(0L),
                                mc.cores=1L,
                                ...){
            
            type <- match.arg(type,eval(formals("plotBias")$type))
            
            .validate(bias=bias)
            
            # dispatch
            p <- .plotBias(bias,type,mc.cores=mc.cores,...)
            
            print(p)
          })

setMethod(f = "plotEqClassQC", 
          signature = c("list","missing"),
          definition = function(
    eqc=list(),
    ...
          ){
            .validate(eqc=eqc)
            
            boxplot(lapply(eq,"[[","size"),
                    las=2,outline=FALSE,
                    main="EQC size distribution per sample")
            
            boxplot(lapply(eqc,"[[","counts"),
                    las=2,outline=FALSE,
                    log="y",
                    main="EQC count distribution per sample")
            
          })

setMethod(f = "plotProportion", 
          signature = c("tbl_df","missing"),
          definition = function(
    tb=tibble(),
    ...
          ){
            .validate(tb=tb)
            
            p <- ggplot(tb %>% pivot_longer(!Transcript,names_to="Sample"),
                        aes(x=Sample,y=value)) + 
              geom_violin() +
              theme(axis.text.x=element_text(angle=90))
            
            plot(p)
            
            invisible(p)
          })

setMethod(f = "plotProportion", 
          signature = c("tbl_df","character"),
          definition = function(
    tb=tibble(),
    tx=character(0L),
    ...
          ){
            .validate(tx=tx,tb=tb)
            
            p <- ggplot(tb %>% 
                     filter(Transcript %in% tx) %>% 
                     pivot_longer(!Transcript,names_to="Sample"),
                   aes(x=Transcript,y=value,color=Sample)) + 
              geom_point() +
              scale_x_discrete("transcript") +
              scale_y_continuous("proportion unique",limits=c(0,1)) +
              theme(axis.text.x=element_text(angle=90))
            
            plot(p)
            
            invisible(p)
          })


#' # Internal functions
".extractBinary" <- function(salmonDir=character(0L),
                             type=character(0L),
                             mc.cores=1L,
                             ...){
  
  obs <- list.files(salmonDir,
                    recursive=TRUE,
                    full.names=TRUE,
                    pattern="obs_gc.gz")
  
  bn <- switch(type,
         "gc" = mapply("-",
                       mclapply(obs,.readGCBinary,mc.cores=mc.cores),
                       mclapply(sub("obs_gc","exp_gc",obs),.readGCBinary,mc.cores=mc.cores),
                       SIMPLIFY=FALSE)
  )
  names(bn) <- basename(dirname(dirname(obs)))
  return(bn)
}

".plotBias" <- function(bias=list(),
                        type=character(0L),
                        mc.cores=1L,
                        simplify=FALSE,
                        ...){
  p <- switch(type,
              "gc" = {
                tb <- do.call(bind_rows,mclapply(names(bias),function(n,bias){
                  b <-  bias[[n]]
                  tibble(sample=n,
                         values=as.vector(b),
                         cond_prob=paste0("cprob",rep(1:nrow(b),ncol(b))),
                         bin=rep(1:ncol(b),each=nrow(b)))
                },bias,mc.cores=mc.cores))
                
                p <- ggplot(tb,
                            aes(y=values,x=bin,group=cond_prob,col=cond_prob)) +
                  ggtitle("Conditional probabilities observed vs. expected GC bias") +
                  scale_y_continuous("log2 fold-change") +
                  scale_x_continuous("bins")
                if ( simplify ){
                  p <- p + geom_smooth()
                } else {
                  p <- p + geom_line() + facet_wrap(~sample)
                }
              })
}

".readGCBinary" <- function(f){
  con <- gzfile(f,"rb")
  lgb <- readBin(con,what="integer",size=4,signed=TRUE)
  nro <- readBin(con,what="integer",size=8,n=1,signed=TRUE)
  nco <- readBin(con,what="integer",size=8,n=1,signed=TRUE)
  mat <- matrix(readBin(con,what="double",size=8,signed=TRUE,n=nro*nco),
                byrow=TRUE,nrow=nro,ncol=nco)
  if ( lgb == 0 ){
    mat <- log2(mat)
  }
  close(con)
  return(mat)
}

".validate" <- function(
    bias=list(),
    dir=character(0L),
    eqc=list(),
    tx=character(0L),
    tb=tibble(),
    type=character(0L)){
  
  if ( length(bias) != 0L ) { 
    stopifnot(length(unique(sapply(bias,nrow)))==1)
    stopifnot(length(unique(sapply(bias,ncol)))==1)
  }
  
  if ( length(dir) != 0L ) { stopifnot(dir.exists(dir)) }
  
  if ( length(eqc) != 0L ) { 
    stopifnot(all(elementNROWS(eqc) == 3))
    stopfinot(all(sapply(lapply(lapply(eqc,names),"%in%",c("size","members","counts")),all)))
  }
  
  if ( nrow(tb) > 0L ) {
    stopifnot("Transcript" %in% colnames(tb))
    stopifnot(ncol(tb) > 1L)
  }
  
  if ( length(tx) != 0L ) { 
    stopifnot(nrow(tb) > 0L) 
    stopifnot("Transcript" %in% colnames(tb))
    stopifnot(all(tx %in% tb$Transcript))
  }
  
  
}
