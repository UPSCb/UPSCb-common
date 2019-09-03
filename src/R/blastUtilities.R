### =======================================================
## libraries
### =======================================================
library(Biostrings)
library(LSD)

### =======================================================
## generics
### =======================================================
setGeneric(name="readBlast",
           def=function(file=character(0),
                        format=character(0),
                        query.fasta=character(0),
                        subject.fasta=character(0),
                        query.name="query",
                        subject.name="subject",
                        plot=TRUE,bestHit=FALSE,
                        ignoreSelf=FALSE,
                        verbose=TRUE){
             standardGeneric("readBlast")
           })

### =======================================================
## methods
### =======================================================
setMethod(f="readBlast",
          signature="character",
          definition=function(
            file=character(0),
            format=c("query.id",
                     "subject.id",
                     "percent.identity",
                     "alignment.length",
                     "mismatches",
                     "gap.opening",
                     "query.start",
                     "query.end",
                     "subject.start",
                     "subject.end",
                     "e.value",
                     "bit.score"),
            query.fasta=character(0),
            subject.fasta=character(0),
            query.name="query",
            subject.name="subject",
            plot=TRUE,bestHit=FALSE,
            ignoreSelf=FALSE,
            verbose=TRUE
            ){
            
            ## check the file size
            if(verbose){
              message("Check file size")
            }
            if(file.info(file)$size == 0){
              warning(paste("The file", file, "has no records"))
              return(NULL)
            }
            
            ## check the format
            if(verbose){ 
              message("Check format")
            }
            firstLine <- read.delim(file=file,nrows=1)
            if(ncol(firstLine)!=length(format)){
              stop("Your format does not match the number of columns in your blast file")
            }
            
            ## define format
            ## 12 entries, we assume blast
            ## more we assume blast+
            fmt <- ifelse(ncol(firstLine)==12,"bm8","blt")
            
            ## read the blast file
            if(verbose){
              message("Reading the file")
            }
            blf <- read.delim(file=file,
                              stringsAsFactors=FALSE,
                              col.names=format)
            
            ## if format is blt, check for extended format
            if(fmt=="blt"){
              if(!all(c("query.length","subject.length") %in% colnames(blf))){
                fmt="bm8"
              }
            }
            
            ## if format is bm8, check for fasta files
            if(fmt=="bm8"){
              ## query cov
              if(length(query.fasta)==1){
                stopifnot(file.exists(query.fasta))
                seq1 <- readDNAStringSet(query.fasta)
                names(seq1) <- sub(" .*","",names(seq1))
                blf$query.length <- width(seq1)[match(blf$query.id,names(seq1))]
              }
              
              ## read the subject
              if(length(subject.fasta)==1){
                stopifnot(file.exists(subject.fasta))
                seq2 <- readDNAStringSet(subject.fasta)
                names(seq2) <- sub(" .*","",names(seq2))
                blf$subject.length <- width(seq1)[match(blf$subject.id,names(seq1))]
              }
            }
            
            ## filter self-self hit
            if(ignoreSelf){
              warning("You have turned ignoreSelf on. This removes any HSPs were query ID == subject ID")
              blf <- blf[blf$query.id != blf$subject.id,]
            }
            
            ## plot the identity
            if(plot){
              plot(density(blf$percent.identity),
                   main=paste(query.name,"-",subject.name," percent identity"))
            }

            ## calculate the percentages
            blf$query.percent <- blf$alignment.length / blf$query.length * 100
            blf$subject.percent <- blf$alignment.length / blf$subject.length * 100
            
            ## plot
            if(plot){
              comparisonplot(blf$query.percent,blf$subject.percent,
                             xlab=paste(query.name,"coverage"),
                             ylab=paste(subject.name,"coverage"),
                             main=paste(query.name,"vs.",subject.name),
                             xlim=range(blf$query.percent),
                             ylim=range(blf$subject.percent))
            }
            
            ## get the cumulative coverage
            ids <- paste(blf$query.id,blf$subject.id,sep="+")
            suids <- sort(unique(ids))
            
            df <- data.frame(query.id=sub("\\+.*","",suids),
                             subject.id=sub(".*\\+","",suids),
                             query.cum.cov=sum(width(IRanges::reduce(split(
                               IRanges(
                                 start=ifelse(blf$query.start>blf$query.end,blf$query.end,blf$query.start),
                                 end=ifelse(blf$query.start<blf$query.end,blf$query.end,blf$query.start)),
                               ids))))/blf$query.length[match(sub("\\+.*","",suids),blf$query.id)],
                             subject.cum.cov=sum(width(IRanges::reduce(split(
                               IRanges(
                                 start=ifelse(blf$subject.start>blf$subject.end,blf$subject.end,blf$subject.start),
                                 end=ifelse(blf$subject.start<blf$subject.end,blf$subject.end,blf$subject.start)),
                               ids))))/blf$subject.length[match(sub(".*\\+","",suids),blf$subject.id)],
                             stringsAsFactors=FALSE)
            
            ## plot again
            if(plot){
              comparisonplot(df$query.cum.cov,df$subject.cum.cov,
                             xlab=paste(query.name,"cumulative coverage"),
                             ylab=paste(subject.name,"cumulative coverage"),
                             main=paste("cumulative coverage",
                                        query.name,"vs.",subject.name),
                             xlim=range(df$query.cum.cov),
                             ylim=range(df$subject.cum.cov))
            }

            ## if best hit
            if(bestHit){
              ord <- order(df$query.cum.cov,decreasing=TRUE)
              sel <- match(unique(df$query.id),df[ord,"query.id"])
              df <- df[ord,][sel,]
              blf <- blf[paste(blf$query.id,blf$subject.id,sep="+") %in% rownames(df),]
            }
            
            ## done
            return(list(blf=blf,df=df))
          })

### =======================================================
## Environment variables
### =======================================================
BM8 <- c("query.id",
         "subject.id",
         "percent.identity",
         "alignment.length",
         "mismatches",
         "gap.opening",
         "query.start",
         "query.end",
         "subject.start",
         "subject.end",
         "e.value",
         "bit.score")
BM8ext <- c(BM8,
            "query.length",
            "subject.length")

### =======================================================
## TODO functionalise if needed
### =======================================================
  
## filter
  
# correct<-FALSE
#   if(!is.null(subsets[[set1]])){
#     sset <- subsets[[set1]]
#     stopifnot(!any(unlist(sset[is.na(match(sset[,2],df$query.id)),2]) %in% df$query.id))
#     if(all(sset[,2]=="")){
#       df <- df[sub("Kmer10_","",df$query.id) %in% sset[,1],]
#       seq1 <- seq1[names(seq1) %in% df$query.id,]
#     }else{
#       sset <- sset[sset[,2] %in% df$query.id,]
#       df[match(sset[,2],df$query.id),"query.id"] <- sset[,1]
#       seq1 <- seq1[!names(seq1) %in% sset[,2],]
#     }    
#     correct<-TRUE
#   }
#   
#   if(!is.null(subsets[[set2]])){
#     sset <- subsets[[set2]]
#     stopifnot(!any(unlist(sset[is.na(match(sset[,2],df$subject.id)),2]) %in% df$subject.id))
#     sset <- sset[sset[,2] %in% df$subject.id,]
#     df[match(sset[,2],df$subject.id),"subject.id"] <- sset[,1]
#     seq2 <- seq2[!names(seq2) %in% sset[,2],]
#     correct<-TRUE
#   }
#   
#   if(correct){
#     df<-df[!duplicated(df[,1:2]),]
#     
#     ## remove the noise
#     df <- df[!(df$query.cum.cov <= 0.3 & df$subject.cum.cov <= 0.3),]
#     comparisonplot(df$query.cum.cov,df$subject.cum.cov,
#                    xlab=paste(set1,"cumulative coverage"),
#                    ylab=paste(set2,"cumulative coverage"),
#                    main=paste("cumulative coverage",set1,"vs.",set2))
#     plot(ecdf(table(df$query.id)),xlim=c(1,20),main=paste(set1,"fragmentation"),xlab="putative number of fragments")
#     q.frag=table(table(df$query.id))/length(unique(df$query.id))
#     plot(ecdf(table(df$subject.id)),xlim=c(1,20),main=paste(set2,"fragmentation"),xlab="putative number of fragments")
#     s.frag=table(table(df$subject.id))/length(unique(df$subject.id))
#   }

## estimate gene number

#   ## find out the isoforms
#   ## careful that we don't drop gene families / paralogs
#   ## if so we will be underestimating.
#   
#   ## select those that are > 0.8 in either direction
#   ##  df8 <- df[df$query.cum.cov >= .8 | df$subject.cum.cov >= 0.8,]
#   
#   ##  df8 <- df8[order(df8$query.cum.cov,decreasing=TRUE),]
#   ##  df8[duplicated(df8$query.id),"subject.id"]
#   
#   ##  df8 <- df8[order(df8$subject.id,decreasing=TRUE),]
#   ##  df8[duplicated(df8$query.id),"subject.id"]
#   
#   ## report
#   cbind(
#     ## 1. the intersect size more than .7 on both axis
#     c(sum(df$query.cum.cov >= 0.7 & df$subject.cum.cov >= 0.7),
#       
#       ## 2. the unique query size
#       length(seq1),
#       
#       ## 3. the unique subject size
#       length(seq2),
#       
#       ## 4. the universe size estimation
#       estimateGeneNumber(length(seq1),length(seq2),sum(df$query.cum.cov >= 0.7 & df$subject.cum.cov >= 0.7))),
#     
#     c(sum(df$query.cum.cov >= 0.7 & df$subject.cum.cov >= 0.7),
#       sum(length(seq1) * q.frag / as.integer(names(q.frag))),
#       sum(length(seq2) * s.frag / as.integer(names(s.frag))),
#       estimateGeneNumber(sum(length(seq1) * q.frag / as.integer(names(q.frag))),
#                          sum(length(seq2) * s.frag / as.integer(names(s.frag))),
#                          sum(df$query.cum.cov >= 0.7 & df$subject.cum.cov >= 0.7)))
#   )
# 
