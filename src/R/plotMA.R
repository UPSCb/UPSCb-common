setMethod(f="plotMA",
          signature="DESeqResults",
          definition=function(object,alpha=0.01,lfc=0.5){
            
            ## lib
            require(LSD)
            
            ## check
            if(!existsFunction("densityPlot")){
              stop("Load the densityPlot function prior to using this function.")
            }
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha & abs(object$log2FoldChange[sel])>=lfc
            
            ## graphic params
            orig.par <- par(no.readonly=TRUE)
            par(mfrow=c(2,1))
            
            ## plots
            par(mar=c(0.1,4.1,4.1,0.1))
            heatscatter(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        add.contour=TRUE,main="MA-plot and MA density estimation",
                        xaxt="n",ylab="log2 FC")
            
            points(log10(object$baseMean[sel][sel2]),
                   object$log2FoldChange[sel][sel2],
                   col="darkred",pch=19,cex=.5)
            
            legend("topright",pch=19,col="darkred","sign. feats.")
            
            par(mar=c(5.1,4.1,0.1,0.1))
            densityPlot(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        grid=250,ncol=30,nlevels=10,
                        main="")
            mtext("log10 mean expression",side=1,line=2)
            title(sub=paste(paste0(sub(".*: ","",elementMetadata(object)$description[2]),":"),
                            sum(sel2),
                            "sign. feats. @",
                            alpha,"FDR and", lfc,
                            "log2FC cutoffs"))
            #mtext("log2 FC",side=2,line=3)
            
            par(orig.par,no.readonly=TRUE)
            invisible(TRUE)
          })
