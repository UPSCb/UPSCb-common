plotSFT <- function(sft,powers=1:35,ymin=NULL,ymax=NULL){
  
  ## check the library
  stopifnot(require(WGCNA))
  
  ## default
  ylim <- range(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  
  ## update ylim
  if(! is.null(ymin)){
    stopifnot(ymin>=-1 & ymin <=1)
    ylim[1] <- ymin
  }
  if(! is.null(ymax)){
    stopifnot(ymax>=-1 & ymax <=1)
    ylim[2] <- ymax
  }
  
  ## do the plot
  plot(0,0,
       xlim=range(sft$fitIndices[,1]),
       ylim=ylim,
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"));
  text(x=sft$fitIndices[,1],
       y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=1,col="red");
  abline(h=0.8 * c(-1,1),col='skyblue');
  abline(h=0.9 * c(-1,1),col='darkolivegreen2')
  
  ## return
  invisible(TRUE)
}
