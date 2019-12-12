#' ---
#' title: "Sample/tissue expression specificity scores"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Rationale
#' This file implements the tissue expression specificity score - which
#' we used to unappropriately call Tau score. The score ranges from 0 to 1;
#' i.e. from ubiquitous to tissue specific and is defined for gene _i_ by:
#' 
#' 1/n-1 * sum(1 - (aij/max(ai in 1:n)))
#' 
#' where _j_ is a given tissue and 
#' _n_ the number of tissues. _aij_ is the average expression of gene _i_ in
#' tissue _j_ and $max(ai in 1:n)$ is the maximal average tissue expression for 
#' gene _i_ accross the _n_ tissues.
#' 
#' The implementation has 2 modes: 'global' and 'local'. The former allows 
#' to calculate the sample specific expression score, _i.e._ is irrelevant of
#' the tissue - which is adequate when the expression cannot be grouped by 
#' tissue, such as is the case with our expression atlases. The second
#' calculates the tissue expression specificity and hence requires several
#' samples per tissue. If the 'local' mode is applied to a dataset such as
#' the expression catalog where every sample represent a different tissue,
#' all expression specificity score will be 0 since every sample is going to
#' maximise its tissue score, hence leading to sum(1-1)/
#' 
#' In the case of missing expression, the expression specificity score is 
#' calculated only for those samples showing expression. _I.e._ _n_ is then
#' the number of samples with expression > 0.
#' 
#' Finally, by defaut the method returns only the specificity scores; _i.e._
#' output="simple", whereas the _complete_ output will report in addition
#' aij, maxj(aij), and n; _i.e_ the number of 
#' samples with expresion > 0 - as well as _peak_, a column that returns the
#' tissue(s) with highest epxression. If multiple, they are comma separated
#' ```{r empty, echo=FALSE, eval=FALSE}
#' ```
#' # Setup
#' Load libraries
suppressPackageStartupMessages(library(matrixStats))

#' # Generic
setGeneric(name="expressionSpecificity",
           def=function(
             exp.mat=matrix(),
             tissues=character(0),
             mode=c("local","global"),
             output=c("simple","complete")
           ){
             standardGeneric("expressionSpecificity")
           })
          
#' # Methods
#' ## matrix-character
setMethod(f="expressionSpecificity",
          signature=c("matrix","character"),
          definition=function(
            exp.mat=matrix(),
            tissues=character(0),
            mode=c("local","global"),
            output=c("simple","complete")
            ){

            # select output
            output <- match.arg(output)
            
            # check
            stopifnot(length(tissues)==ncol(exp.mat))
            
            # turn no expression into NA
            exp.mat[exp.mat==0]<-NA
            
            # switch based on the mode
            # to calculate the numerator and denominator
            frac <- switch(
              match.arg(mode),
              local = {
                message("Calculating tissue specificity")
                list(
                  aij = do.call(cbind,
                                 lapply(split.data.frame(t(exp.mat),
                                                         tissues,
                                                         drop=FALSE),
                                        colMeans,na.rm=TRUE)),
                  maxn = rowMaxs(do.call(cbind,
                                  lapply(split.data.frame(t(exp.mat),
                                                          tissues,
                                                          drop=FALSE),
                                         colMeans,na.rm=TRUE)),na.rm=TRUE))
              },
              global = {
                message("Calculating sample specificity")
                list(aij = exp.mat,
                     maxn = rowMaxs(exp.mat,na.rm=TRUE))
              }
            )
            
            # localSums
            vals <- 1 - (frac$aij/frac$maxn)
            
            # account for no expression
            svals <- rowSums(vals,na.rm=TRUE)
            n <- rowSums(frac$aij != 0,na.rm=TRUE)
            
            # get the score
            scores <- svals / (n -1)
            
            # since n-1 can be 0 and that would mean that sample
            # are only expressed in one sample, we just convert NaN to 1
            scores[is.nan(scores)] <- 1
            
            ## and turn the genes without expression to NA
            scores[n==0] <- NA
            
            # expressionSpecificity score
            if(output=="simple"){
              return(scores)
            } else {
              return(data.frame(
                score=scores,
                aij=frac$aij,
                maxn=frac$maxn,
                n=n,
                peak=sapply(apply(frac$aij == frac$maxn,1,subset,x=colnames(frac$aij)),paste,collapse=",")
              ))
            }
            })


#' ## data-frame character
#' It simply dispatch to matrix-character
setMethod(f="expressionSpecificity",
          signature=c("data.frame","character"),
          definition=function(
            exp.mat=data.frame(),
            tissues=character(0),
            mode=c("local","global"),
            output=c("simple","complete")
            ){
            expressionSpecificity(
              exp.mat=as.matrix(exp.mat),
              tissues=tissues,
              mode=mode,
              output=output
              )
          })

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
