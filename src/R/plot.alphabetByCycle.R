setGeneric(name="plotAlphabetByCycle",
          def=function(obj=DNAStringSet(),...){
            standardGeneric("plotAlphabetByCycle")
          })

setMethod(f="plotAlphabetByCycle",
          signature="DNAStringSet",
          definition=function(obj,smooth=FALSE,...){
            
            ## some failsafe checks
            stopifnot(require(RColorBrewer))
            stopifnot(require(ShortRead))
            stopifnot(length(obj)>0)
            stopifnot(length(unique(width(obj)))==1)
            
            ## abc
            abc <- alphabetByCycle(obj)[c("A","C","G","T"),]
            
            if(smooth){
              abc <- t(sapply(lapply(apply(abc,1,Rle),
                                     runmean,k=9,endrule="constant"),
                              as.integer))
            }
            
            ## plot it
            rng <- range(abc)
            pal <- brewer.pal(4,"Dark2")
            plot(0,0,type="n",ylim=rng,xlim=c(0,ncol(abc)),
                 ylab="occurence",xlab="bp/cycle")
            sapply(1:4,function(i,abc,pal,...){
              lines(abc[i,],col=pal[i],...)},abc,pal,...)
            legend("top",lty=1,col=pal,c("A","C","G","T"),...)
          })

setMethod(f="plotAlphabetByCycle",
          signature="AAStringSet",
          definition=function(obj,log="",
                              legend.x="top",
                              IUPAC=FALSE,...){
            
            ## fail safe checks
            stopifnot(require(RColorBrewer))
            stopifnot(require(ShortRead))
            
            ## there's a lot hardcoded here
            ## no check for the aa            
            ## no check for the color            
            if(IUPAC){
              alph <- names(AMINO_ACID_CODE)
            } else {
              alph <- names(AMINO_ACID_CODE)[1:20]
            }
            
            aa.mat <- do.call(rbind,strsplit(as.character(obj),""))
            abc <- apply(aa.mat,2,function(co,alph){
              table(factor(co,levels=alph))},alph)
            
            ## add a pseudo count to avoid
            ## infinite log issue
            if(log=="y"){
              abc <- abc+1
            }
            
            ## plot it
            rng <- range(abc)
            pal <- rep(brewer.pal(8,"Dark2"),3)
            ltys <- rep(1:3,each=8)
            plot(0,0,type="n",ylim=rng,xlim=c(0,ncol(abc)),
                 ylab="occurence",xlab="aa/cycle",log=log,...)
            
            sapply(1:nrow(abc),function(i,abc,pal,ltys,...){
              lines(abc[i,],col=pal[i],lty=ltys[i],...)},abc,pal,ltys,...)
            legend(legend.x,col=pal,lty=ltys,rownames(abc),...)
          })
