pcaPlotTsd = function(x,comp1=1,comp2=2,screeplot=FALSE, adj.lim=c(0.001,0.1),
                     treatment=treatment,sample.ids=sample.ids,context,
                     scale=TRUE,center=TRUE,obj.return=FALSE){

  x.pr = prcomp(t(x),scale.=scale,center=center)
  x.pr.pv = round(x.pr$sdev^2/sum(x.pr$sdev^2),2)*100


  if (screeplot){
    i=5;screeplot(x.pr, type="barplot",
                  main=paste(context,"methylation PCA Screeplot"),
                  col = rainbow(i)[i])
  }
  else{
    #loads = loadings(x.pr)
    #loads = x.pr$rotation
    treatment=treatment
    sample.ids=sample.ids
    my.cols=rainbow(length(unique(treatment)), start=1, end=0.6)
    pc1=x.pr$x[,comp1]
    pc2=x.pr$x[,comp2]

    plot(pc1,pc2, main = paste(context,"methylation PCA Analysis"),
         col=my.cols[treatment+1],
         xlab=paste(x.pr.pv[1],"%"), ylab=paste(x.pr.pv[2],"%"))
    text(pc1, pc2,labels=sample.ids,adj=c(-0.4,0.3), col=my.cols[treatment+1])
  }
  if(obj.return){  return((x.pr))}
}


newPCASamples <- function (.Object, screeplot = FALSE, adj.lim = c(4e-04, 0.1),
          scale = TRUE, center = TRUE, comp = c(1, 2), transpose = TRUE,
          sd.filter = TRUE, sd.threshold = 0.5, filterByQuantile = TRUE,
          obj.return = FALSE)
{
  require(matrixStats)
  mat = getData(.Object)
  mat = mat[rowSums(is.na(mat)) == 0, ]
  meth.mat = mat[, .Object@numCs.index]/(mat[, .Object@numCs.index] +
                                           mat[, .Object@numTs.index])
  names(meth.mat) = .Object@sample.ids
  if (sd.filter) {
    if (filterByQuantile) {
      sds = rowSds(as.matrix(meth.mat))
      cutoff = quantile(sds, sd.threshold)
      meth.mat = meth.mat[sds > cutoff, ]
    }
    else {
      meth.mat = meth.mat[rowSds(as.matrix(meth.mat)) >
                            sd.threshold, ]
    }
  }
  if (transpose) {
    pcaPlotTsd(meth.mat, comp1 = comp[1], comp2 = comp[2],
              screeplot = screeplot, adj.lim = adj.lim, treatment = .Object@treatment,
              sample.ids = .Object@sample.ids, context = .Object@context,
              scale = scale, center = center, obj.return = obj.return)
  }
  else {
    .pcaPlot(meth.mat, comp1 = comp[1], comp2 = comp[2],
             screeplot = screeplot, adj.lim = adj.lim, treatment = .Object@treatment,
             sample.ids = .Object@sample.ids, context = .Object@context,
             scale = scale, center = center, obj.return = obj.return)
  }
}