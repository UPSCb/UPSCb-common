#' ---
#' title: "Feature selection"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions to 
#' 
#' 1. select feature
#' provided they are expressed at an `exp` cutoff (exp) in at least `nrep` 
#' replicates of any `condition`.
#' 
#' 2. plot 1. over a range of cutoffs
#' 
#' 3. generate samples summary of how many genes are expressed using the
#' given cutoff
#' 
#' 4. a function to plot 3.
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/featureSelection.R")
#' ```
#' 
#' In the following, you will first find the S4 generic of the function
#'  and then its implementation, as well as some description of its arguments.
#' 
#' # Generic
setGeneric(name="featureSelect",
           def=function(counts=matrix(),
                        conditions=factor(),
                        exp=1L,
                        nrep=2L,
                        scale=FALSE
           ){
             standardGeneric("featureSelect")
           })

setGeneric(name="rangeFeatureSelect",
           def=function(counts=matrix(),
                        conditions=factor(),
                        nrep=2L,
                        plot=TRUE,
                        scale=FALSE
           ){
             standardGeneric("rangeFeatureSelect")
           })

setGeneric(name="featureSelectProp",
           def=function(counts=matrix(),
                        conditions=factor(),
                        exp=1
           ){
               standardGeneric("featureSelectProp")
           })

setGeneric(name="samplesSummary",
           def=function(counts=matrix(),
                        conditions=factor(),
                        exp=1L,
                        nrep=2L
           ){
               standardGeneric("samplesSummary")
           })

setGeneric(name="rangeSamplesSummary",
           def=function(counts=matrix(),
                        conditions=factor(),
                        nrep=2L,
                        plot=TRUE
           ){
               standardGeneric("rangeSamplesSummary")
           })


#' # Implementation of featureSelect
#' 
#' ## Arguments
#' The 'featureSelect' function expect the following arguments:
#' 
#' 1. `counts` the matrix of expression counts (normalised or not)
#' 
#' 2. `conditions` a factor that describes the different condition in the 
#' experiment
#'
#' It has 2 optional arguments
#' 
#' 1. `exp` the minimal value of expression to be used as a cutoff
#' 
#' 2. `nrep` the minimal number of replicates that need to have an expression
#' over the `exp` cutoff in any given condition
#'
#' ## Value
#' 
#' The function returns a boolean vector where TRUE indicates that the corresponding
#' feature passed the selected thresholds. The features are ordered as the 
#' rows of the input matrix
#'
#' # Implementation of featureSelectProp
#'
#' ## Arguments
#' The 'featureSelectProp' function expect the following arguments:
#' 
#' 1. `counts` the matrix of expression counts (normalised or not)
#' 
#' 2. `conditions` a factor that describes the different condition in the 
#' experiment
#'
#' It has 1 optional argument
#' 
#' 1. `exp` the proportion of reads within one set of replicates used as a
#' cutoff
#' 
#' ## Value
#' 
#' The function returns a boolean vector where TRUE indicates that the corresponding
#' feature passed the selected thresholds. The features are ordered as the 
#' rows of the input matrix

setMethod(f = "featureSelect", 
          signature = c("matrix","factor"),
          definition = function(
            counts=matrix(),
            conditions=factor(),
            exp=1L,
            nrep=2L,
            scale=FALSE
          ){
            
            ## validation
            .validate(counts,conditions)
             
            ## compute
            rowSums(.count(counts=counts,
                            conditions=conditions,
                            exp=exp,nrep=nrep,scale=scale),na.rm=TRUE) >= 1
                        
          })

setMethod(f = "rangeFeatureSelect", 
          signature = c("matrix","factor"),
          definition = function(
            counts=matrix(),
            conditions=factor(),
            nrep=2L,
            plot=TRUE,
            scale=FALSE){
              
            ## validation
              .validate(counts,conditions)
            
            ## compute
            rg <- switch(as.character(scale),
                "FALSE" = 0:floor(max(counts)),
                "TRUE" = seq(0,2,0.1)
            )
            res <- lapply(rg,function(i){
              featureSelect(counts=counts,
                            conditions=conditions,
                            exp=i,nrep=nrep,scale=scale)
            })
            
            ## plot
            if(plot){
              plot(rg,sapply(res,sum),
                   main="Number of selected genes by cutoff values",
                   type="b",ylab="Number of genes",xlab="cutoff")
                #,xaxt="n")
              #axis(1,rg+1,labels=rg)
              
              plot(rg,sapply(res,sum),log="y",
                   main="Number of selected genes by cutoff values",
                   type="b",ylab="Number of genes",xlab="cutoff")
              #,xaxt="n")
              #axis(1,rg+1,labels=rg)
            }
            
            ## return
            return(res)
          })

setMethod(f = "featureSelectProp", 
          signature = c("matrix","factor"),
          definition = function(
              counts=matrix(),
              conditions=factor(),
              exp=1
          ){
              
              ## validation
              .validate(counts,conditions)
              
              ## compute
              rowSums(sapply(lapply(
                  split.data.frame(t(sweep(counts, 2, colSums(counts), "/")),conditions)
                  ,colSums), ">=", exp)) >= 1
              
          })

setMethod(f = "samplesSummary", 
          signature = c("matrix","factor"),
          definition = function(
              counts=matrix(),
              conditions=factor(),
              exp=1L,
              nrep=2L
          ){
              
              ## validation
              .validate(counts,conditions)
              
              ## compute
              colSums(.count(counts=counts,
                              conditions=conditions,
                              exp=exp,nrep=nrep))
              
          })

setMethod(f = "rangeSamplesSummary", 
          signature = c("matrix","factor"),
          definition = function(
              counts=matrix(),
              conditions=factor(),
              nrep=2L,
              plot=TRUE){
              
              ## validation
              .validate(counts,conditions)
              
              ## compute
              mx <- floor(max(counts))
              res <- sapply(0:mx,function(i){
                  samplesSummary(counts=counts,
                                conditions=conditions,
                                exp=i,nrep=nrep)
              })
              
              ## plot
              if(plot){
                  
                  heatmap.2(log10(res+1),dendrogram="row",
                            Colv=FALSE,
                            col=colorRampPalette(c("blue","white","red"))(100),
                            main="log10 # of genes",xlab="vst expression cutoff")
                  
                  meds <- colMedians(res) + 1
                  barplot(t(t(res) / meds) - 1,
                          beside=TRUE,col=rainbow(nlevels(conditions)),
                          main="# gene expressed / median # gene expressed")
                  legend("topleft",legend=levels(conditions),fill=rainbow(nlevels(conditions)))
              }
              
              ## return
              return(res)
          })

# internal function for summarising
.count <- function(counts=matrix(),
                   conditions=factor(),
                   exp=1L,
                   nrep=2L,scale=FALSE){
    counts <- switch(as.character(scale),
        "FALSE" = t(counts),
        "TRUE" = scale(t(counts))
    )
    sapply(lapply(
        split.data.frame(counts >= exp,conditions)
        ,colSums), ">=", nrep)
}

.validate <- function(counts,conditions){
    stopifnot(all(dim(counts) > 1))
    stopifnot(nlevels(conditions)>=2)
}
