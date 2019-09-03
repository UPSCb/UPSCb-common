setGeneric(name="readPfam",def=function(file=character(0)){
  standardGeneric("readPfam")
})

setMethod(f="readPfam",
          signature="character",
          definition=function(file){
            
            ## read the file
            PFAM <- read.delim(gzfile(file),
                       stringsAsFactors=FALSE,header=FALSE, 
                       col.names=c("geneID","PFAMID","Desc"))
            
            ## change relationship from one to many to one to one
            ## i.e. gene ids get duplicated
            PFAM$Desc[grep("PF13720",PFAM$PFAMID)] <- sub("; Domain 2","",PFAM$Desc[grep("PF13720",PFAM$PFAMID)])
            
            ## Change all description to lower case
            PFAM$Desc <- tolower(PFAM$Desc)
            
            ## return
            data.frame(geneID=rep(PFAM$geneID,sapply(strsplit(PFAM$PFAMID,";"),length)),
                       PFAMID=unlist(strsplit(PFAM$PFAMID,";")),
                       Desc=unlist(strsplit(PFAM$Desc,";")))
          })

setGeneric(name="extractPfamId",def=function(obj,set){
  standardGeneric("extractPfamId")
})

setMethod(f="extractPfamId",
          signature=c("data.frame","character"),
          definition=function(obj,set){
  
            ## a series of check
            stopifnot(all(c("geneID","PFAMID") %in% colnames(obj)))
            #stopifnot(all(set %in% obj$geneID))
            
            ## get the list
            obj[obj$geneID %in% set,"PFAMID"]
          })

setMethod(f="extractPfamId",
          signature=c("data.frame","factor"),
          definition=function(obj,set){
            extractPfamId(obj,as.character(set))
          })

setGeneric(name="pfamEnrichment",def=function(set,population,alpha=0.05,plot=TRUE){
  standardGeneric("pfamEnrichment")
})

setMethod(f="pfamEnrichment",
          signature=c("factor","factor"),
          definition=function(set,population,alpha=0.05,plot=TRUE){

            ## a series of check
            stopifnot(all(grepl("PF", population)))
            stopifnot(all(grepl("PF", set)))
            stopifnot(all(set %in% population))
            
            ## tabulate the population and set
            pop <- table(as.character(population))
            tset <- table(as.character(set))
            
            ## calculate p.vals
            p.vals <- sapply(as.character(set),function(nam,pop,samp){
              fisher.test(matrix(c(samp[nam],sum(samp) - samp[nam],
                                   pop[nam],sum(pop) - pop[nam]),
                                 ncol=2,byrow=TRUE))$p.value
            },pop,tset)
            names(p.vals) <- set
            
            ## plot
            if(plot){
              hist(unlist(p.vals),breaks=seq(0,1,.01),
                   main="PFAM enrichment p.value",
                   xlab="p. value")
            }
            
            ## calculate adj. p. vals
            adj.p <- p.adjust(p.vals,method="BH")
            if(plot){
              hist(unlist(adj.p),breaks=seq(0,1,.01),
                   main="PFAM enrichment adj. p.value",
                   xlab="adjusted p. value")
            }
            
            ## return the enriched PFAM
            as.character(set[adj.p<=alpha])
            
          })

setGeneric(name="pfamCloud",def=function(obj,set,plot=TRUE,
                                         random.color=TRUE,
                                         color=TRUE,...){
  standardGeneric("pfamCloud")
})

setMethod(f="pfamCloud",
          signature=c("data.frame","character"),
          definition=function(obj,set,
                              plot=TRUE,
                              random.color=TRUE,
                              color=TRUE,...){
            
            ## load lib
            stopifnot(require(wordcloud))
            
            if(color){
              pal=brewer.pal(8,"Dark2")
            } else {
              pal="black"
            }

            ## some checks
            stopifnot(all(c("PFAMID","Desc") %in% colnames(obj)))
            stopifnot(all(set %in% obj$PFAMID))
            
            ## wordcloud
            desc <- obj[match(set,obj$PFAMID),"Desc"]
            
            tab <- sort(table(as.character(desc)),decreasing=TRUE)
            
            
            if(plot){
              wordcloud(names(tab),
                        freq=tab,
                        random.color=random.color,
                        colors=pal,...)
            }
            
            ## return
            invisible(desc)
            })
