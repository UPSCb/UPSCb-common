"WgcnaClusterPlot" <- function(dat,xlabels=NULL,...){
  
  ## check
  stopifnot(require(LSD))
  
  ## matrix
  stopifnot(is.matrix(dat))
  
  ## labels
  if(is.null(xlabels)){
    message("No labels were provided, defining them")
    #     if(ncol(dat)>nrow(dat)){
    #       message("Setting the row names as labels")
    #       xlabels <- rownames(dat)
    #     } else {
    xlabels <- colnames(dat)
    message("Setting the col names as labels")
    #     }
  }
  
  ## plot
  clusterplot(dat,
              colpal=colorRampPalette(c("dodgerblue3",
                                        "lightcyan3",
                                        "lightgrey"))(9),
              quartiles.col=c("black","darkgrey","darkgrey"),
              xlabels=xlabels,ylab="vst expression",xlab="",...)
  
  ## return
  invisible(TRUE)
}