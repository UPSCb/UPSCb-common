#' ---
#' title: "Multi density plot"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

stopifnot(suppressPackageStartupMessages(require(RColorBrewer)))

setGeneric(name="plot.multidensity",
           def=function(x,
                        xlab="x",col=brewer.pal(8,"Dark2"),
                        legend.x="top",xlim=NULL,ylim=NULL,
                        lty=1,legend.cex=1,legend=NULL,
                        legend.col=NULL,legend.lwd=1,...){
  standardGeneric("plot.multidensity")
})

setMethod(f="plot.multidensity",
          signature=c("list"),
          definition=function(x,
                              xlab="x",col=brewer.pal(8,"Dark2"),
                              legend.x="top",xlim=NULL,ylim=NULL,
                              lty=1,legend.cex=1,legend=NULL,
                              legend.col=NULL,legend.lwd=1,...){
            
            ## densities
            dens <- lapply(x,density)
            .plot.multidensity(dens,xlab,col,legend.x,xlim,ylim,
                               lty,legend.cex,legend,legend.col,legend.lwd,...)
          })

setMethod(f="plot.multidensity",
          signature=c("matrix"),
          definition=function(x,
                              xlab="x",col=brewer.pal(8,"Dark2"),
                              legend.x="top",xlim=NULL,ylim=NULL,
                              lty=1,legend.cex=1,legend=NULL,
                              legend.col=NULL,legend.lwd=1,...){
            
            ## densities
            dens <- apply(x,2,density)
            .plot.multidensity(dens,xlab,col,legend.x,xlim,ylim,
                               lty,legend.cex,legend,legend.col,legend.lwd,...)
          })

".plot.multidensity" <- function(d, xlab="x",col=brewer.pal(8,"Dark2"),
                                 legend.x="top",xlim=NULL,ylim=NULL,
                                 lty=1,legend.cex=1,legend=NULL,
                                 legend.col=NULL,legend.lwd=1,...){
  
  if(is.null(xlim)){
    xlim <- range(sapply(d,"[[","x"))
  }
  
  if(is.null(ylim)){
    ylim <- range(sapply(d,"[[","y"))
  }
  
  ## lty
  if(length(lty)==1){
    lty <- rep(lty,length(d))
  }
  
  ## plot
  plot(0,0,xlim=xlim,ylim=ylim,ylab="density",type="n",xlab=xlab,...)
  
  ## lines
  sapply(1:length(d),function(i,d,col,...){lines(d[[i]],col=col[i],lty=lty[i],...)},d,col,...)
  
  ## legend
  if(!is.null(legend) | !is.null(names(d))){
    if(is.null(legend)){
      leg <- names(d)
    } else {
      leg <- legend
    }
    if(!is.null(legend.col)){
      lcol <- legend.col
    } else {
      lcol <- col[1:length(d)]
    }
    legend(legend.x,col=lcol,bty="n",
           legend=leg,lty=lty,cex=legend.cex,
           lwd=legend.lwd)
  }
}
