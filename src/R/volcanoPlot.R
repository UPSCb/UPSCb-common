require(limma)

setGeneric(name="volcanoPlot",def=function(object,alpha=0.01,lfc=0.5,coef=1){
  standardGeneric("volcanoPlot")
})

setMethod(f="volcanoPlot",
          signature="DataFrame",
          definition=function(object,alpha=0.01,lfc=0.5){

            ## lib
            require(LSD)

            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha & 
              abs(object$log2FoldChange[sel]) >= lfc

            ## plot
            heatscatter(object$log2FoldChange[sel],
                        -log10(object$padj[sel]),
                        main="Volcano",xlab=elementMetadata(object)$description[2],
                        ylab=paste("- log(10)",elementMetadata(object)$description[6]))

            ## legend
            legend("top",bty="n",paste("cutoff @",alpha,"and",lfc),lty=2,lwd=2,col="gray")

            ## points
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="lightblue",pch=19)
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="dodgerblue3",pch=19,cex=0.5)

            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,lwd=2,col="gray")
            if(lfc != 0) {
              abline(v=lfc*c(-1,1),lty=2,lwd=2,col="gray")
            }
            })


setMethod(f="volcanoPlot",
          signature="MArrayLM",
          definition=function(object,alpha=0.01,lfc=0.5,coef=1){

            ## lib
            require(LSD)
            require(limma)
            #get topTable
            tt <- limma::toptable(object, coef = coef,
                                  number = length(object$p.value))

            ## selectors
            sel <- tt$adj.P.Val <= alpha

            ## plot
            heatscatter(tt$logFC,
                        -log10(tt$adj.P.Val),
                        main="Volcano",xlab="log2 Fold Change",
                        ylab=paste("- log(10)"," Adjusted p-value"))

            ## legend
            legend("top",bty="n",paste("cutoff @",alpha,"and",lfc),lty=2,lwd=2,col="gray")

            ## points
            points(tt$logFC[sel],-log10(tt$adj.P.Val[sel]),col="lightblue",pch=19)
            points(tt$logFC[sel],-log10(tt$adj.P.Val[sel]),col="dodgerblue3",pch=19,cex=0.5)

            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,lwd=2,col="gray")
          })

setMethod(f="volcanoPlot",
          signature="data.frame",
          definition=function(object,alpha=0.01,lfc=0.5){
            
            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$adj.P.Val)
            sel2 <- object$adj.P.Val[sel]<=alpha & 
              abs(object$logFC[sel]) >= lfc
            
            ## plot
            heatscatter(object$logFC[sel],
                        -log10(object$adj.P.Val[sel]),
                        main="Volcano plot",xlab="log2 Fold Change",
                        ylab="- log(10) BH adjusted p-value")
            
            ## legend
            legend("top",bty="n",paste("cutoff @",alpha,"and",lfc),lty=2,lwd=2,col="gray")
            
            ## points
            points(object$logFC[sel][sel2],-log10(object$adj.P.Val[sel][sel2]),col="lightblue",pch=19)
            points(object$logFC[sel][sel2],-log10(object$adj.P.Val[sel][sel2]),col="dodgerblue3",pch=19,cex=0.5)
            
            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,lwd=2,col="gray")
            if(lfc != 0) {
              abline(v=lfc*c(-1,1),lty=2,lwd=2,col="gray")
            }
          })

