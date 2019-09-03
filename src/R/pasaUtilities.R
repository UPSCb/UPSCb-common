#' ---
#' title: "PASA utilities"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions that help manipulate
#' various PASA result files.
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/pasaUtilities.R")
#' ```
#' 
#' In the following, you will first find the S4 generics of the different 
#' functions and then their implementation, as well as some description of
#' their functionality, arguments, etc.
#' 
#' # Libraries
#' 
#' The following libraries are necessary for the correct execution of the utility.
#' As this file is meant to be _sourced_, we _require_ the presence of the 
#' libraries.
stopifnot(suppressPackageStartupMessages(require(genomeIntervals)))
stopifnot(suppressPackageStartupMessages(require(IRanges)))

#' # Other helper
source("~/Git/UPSCb/src/R/gff3Utilities.R")

#' 
#' # Generics
setGeneric(name="checkAndUpdatePasaGff3IDs",
           def=function(gff3=character(0),
                        output=character(0),
                        ordered.types=character(0),
                        src=NA,mc.cores=8,
                        id.type=c("incremental","MakerP"),
                        rm.dup=TRUE,fix=TRUE
                        ){
             standardGeneric("checkAndUpdatePasaGff3IDs")
           })

setGeneric(name="validateProcessedGff3",
           def=function(
             gff3=character(0),
             types=character(0),
             hierarchy=matrix()
           ){
             standardGeneric("validateProcessedGff3")
           })

#' # Implementation
#' ## checkAndUpdatePasaGff3IDs
#' This function takes up to 4 arguments:
#' 
#' 1. gff3: the input file
#' 
#' 2. output: the output file
#' 
#' 3. ordered.types: the values of the gff3 type column in their desired order
#' for sorting. I.e. gene should come before mRNA, etc.
#' 
#' 4. src: used to fill the output gff3 'source' column for changed genes
#' 
#' 5. mc.cores: use that many cpu for parallel processing
#' 
#' 6. id.type: the type of ID to generate. One of _incremental_ or _MakerP_, 
#' where for the former, the newID is simply the existing highest ID +1, and for
#' the latter, the newID is calculated based on existing IDs, e.g. if a newID
#' needs created between gene 00010 and 00020, then 00015 is selected. _incremental_
#' is the safest and the default.
#' 
#' This function write the gff3 file with corrected IDs into the file provided as
#' 'output' and returns a data.frame that contains the gene ID lift-over. All the
#' other IDs (mRNA, exon, etc.) have been re-assigned using standard conventions
#' for the modified gene IDs (and these only).
setMethod(f = "checkAndUpdatePasaGff3IDs", 
          signature = c("character","character"),
          definition = function(
            gff3=character(0),
            output=character(0),
            ordered.types=c("gene","mRNA","exon","five_prime_UTR","CDS","three_prime_UTR"),
            src=NA, mc.cores=8,
            id.type=c("incremental","MakerP"),
            rm.dup=TRUE,fix=TRUE
          ){
            
            totalSteps=8
            ## ------------
            ## INITIALISE
            ## ------------
            message(sprintf("Initializing 1/%s",totalSteps))
            # some checks
            stopifnot(length(gff3)>0)
            stopifnot(file.exists(gff3))
            stopifnot(length(output)>0)
            
            id.type <- match.arg(id.type)
            
            # then proceed
            gff <- readGff3(gff3,quiet=TRUE)
            message(sprintf("Read %s entries",nrow(gff)))
            
            # remove absolute duplicates
            dup.sel <- duplicated(cbind(gff,annotation(gff)))
            message(sprintf("Removing %s exact duplicates",sum(dup.sel)))
            gff <- gff[!dup.sel,]
            
            # reorder the gff
            typs <- gff$type
            for(l in rev(ordered.types)){typs <- relevel(typs,l)}
            gff$type <- typs
            
            # by scaffold, coordinates, strand and type
            # strand is done by selecting the star or the end for the second
            # item in the selector below
            # gff <- gff[order(seqnames(gff),
            #                 #ifelse(gff$strand=="+",
            #                        gff[,1],
            #                #        gff[,2]),
            #                 gff$type
            #                 ),]
            
            # extract all IDs
            IDs <- getGffAttribute(gff,"ID")
            
            # extract all Parents
            Parents <- getGffAttribute(gff,"Parent")
            
            # select the genes and mRNA
            sel <- gff$type == "gene"
            mRNA.sel <- gff$type == "mRNA"
            
            # if incremental, find the max ID
            # from non novel genes
            nov.sel <- grepl("novel_gene",IDs[sel])
            maxID <- ifelse(id.type=="incremental",
                            max(suppressWarnings(
                              as.integer(
                                substr(IDs[sel][!nov.sel],
                                       nchar(IDs[sel][!nov.sel])-4,
                                       nchar(IDs[sel][!nov.sel])))),na.rm = TRUE),
                            0)
            
            ## ------------
            ## CURATION
            ## ------------
            message(sprintf("Curating mRNA IDs 2/%s",totalSteps))
            # make sure that gene duplication is not just 2 transcripts
            dupIDs <- IDs[mRNA.sel][duplicated(IDs[mRNA.sel])]
            toCurate <- names(table(dupIDs)[table(dupIDs) != 1])
            
            if(length(toCurate)>0){
            
              # check for mRNA overlap
              grngList <- split(as(gff[mRNA.sel][IDs[mRNA.sel] %in% toCurate],"GRanges"),
                                IDs[mRNA.sel][IDs[mRNA.sel] %in% toCurate])
              
              # update mRNA ID
              for(mID in names(elementNROWS(reduce(grngList)))[elementNROWS(reduce(grngList))==1]){
                # rename
                all.pos <- which(Parents==mID & !is.na(Parents))
                i=0
                for (m.pos in rev(which(IDs==mID)[-1])){
                  p.sel <- all.pos[all.pos > m.pos]
                  Parents[p.sel] <- paste0(Parents[p.sel],i)
                  IDs[m.pos] <- paste0(IDs[m.pos],i)
                  i=+1
                  all.pos <- all.pos[all.pos < m.pos]
                }
              }
            
              ## check for mRNA with different parents
              d.sel <- which(IDs[mRNA.sel] %in% dupIDs)
              tab <- sapply(split(Parents[mRNA.sel][d.sel],
                                  IDs[mRNA.sel][d.sel]),unique)
              toCurate <- tab[elementNROWS(tab) > 1]
              
              # update the IDs with the scaffold
              p.sel <- Parents %in% names(toCurate) & ! is.na(Parents)
              Parents[p.sel] <- paste(seqnames(gff)[p.sel],Parents[p.sel],sep="_")
              m.sel <- IDs %in% names(toCurate)
              IDs[m.sel] <- paste(seqnames(gff)[m.sel],IDs[m.sel],sep="_")
            }
                
            # dups on the same scaffold
            dupIDs <- IDs[mRNA.sel][duplicated(IDs[mRNA.sel])]
            d.sel <- which(IDs[mRNA.sel] %in% dupIDs)
            tab <- sapply(split(Parents[mRNA.sel][d.sel],
                                IDs[mRNA.sel][d.sel]),unique)
            toCurate <- tab[elementNROWS(tab) > 1]
            
            if(length(toCurate)>0){
              for ( n in names(toCurate)){
                g.sel <- which(IDs==n)
                p.sel <- which(Parents==n & !is.na(Parents))
                
                Parents[p.sel] <- 
                  paste(seqnames(gff)[g.sel][1],
                        Parents[g.sel][rowSums(sapply(g.sel,"<",p.sel))],
                        "mRNA","1",sep="_")
                IDs[g.sel] <- 
                  paste(seqnames(gff)[g.sel][1],
                        Parents[g.sel],
                        "mRNA","1",sep="_")
              }
            }
            
            # create a dataframe to return gene info
            df <- data.frame(oldID=IDs[sel],
                             newID=as.character(IDs[sel]),
                             status="stable",
                             scaffold=seqnames(gff)[sel],
                             stringsAsFactors = FALSE)
            
            ## ------------
            ## GENE SPLIT
            ## ------------
            ## CAVEAT here we hope that there is no novel gene around
            ## as we do not take this issue in consideration. It might 
            ## mean to immediately extend the gff, and recompute/extend
            ## the various selectors and the ID lists
            message(sprintf("Considering gene split 3/%s",totalSteps))
            # check for duplicated mRNA IDs, i.e. they originate from split genes
            split.gff = df.split = se.rm = NULL
            if(any(duplicated(IDs[mRNA.sel]))){
              
              dupIDs <- IDs[mRNA.sel][duplicated(IDs[mRNA.sel])]
              
              df[df$newID %in% unique(Parents[IDs %in% dupIDs]),"status"] <- "split"
              
              # get the unique dups
              for(pos in match(unique(dupIDs),IDs[mRNA.sel])){
              
                id <- Parents[mRNA.sel][pos]
              
                message(sprintf("Fixing putative split gene at position %s with the ID %s",pos,id))

                # get the entries spanned by that gene
                spanPos <- .getPosition("gene",IDs,Parents,id)

                # get the gene ID
                prev <- df[df$newID == id,]
                
                # stop if not incremental
                if(id.type != "incremental"){
                  stop("DEFUNCT - only incremental IDs are supported - select id.type='incremental'")
                }
                # create a gff3
                tmp.gff <- gff[spanPos,]
                
                # find the number of mRNA and gene
                m.sel <- tmp.gff$type=="mRNA"
                mRNA.num <- sum(m.sel)
                oldmID <- getGffAttribute(tmp.gff[m.sel],"ID")[,1]
                m.grng <- as(tmp.gff[m.sel,],"GRanges")
                g.grng <- reduce(m.grng)
                gene.num <- length(g.grng)
                
                # if we have only one gene; it's simply a duplicated mRNA.ID
                if(gene.num==1){
                  message("Skipping as this is not a split.")
                  df[df$oldID==prev$oldID,"status"] <- "stable"
                  next
                }
                
                # sort it by coordinates, independently of the strand
                # tmp.gff <- tmp.gff[order(tmp.gff[,1],tmp.gff$type),]
                
                # create new genes
                tmp.gff <- c(tmp.gff[rep(which(tmp.gff$type=="gene"),gene.num-1),],tmp.gff)
                
                # update the gene coordinates
                tmp.gff[1:gene.num,1] <- start(g.grng)
                tmp.gff[1:gene.num,2] <- end(g.grng)
                
                # the gene attributes
                maxIDs <- (maxID+1):(maxID + gene.num)
                newIDs <- sprintf("%sg%05d",seqnames(tmp.gff[1:gene.num,]),maxIDs)
                tmp.gff[1:gene.num,]$gffAttributes <- paste("ID=",newIDs,";Name=",newIDs,sep="")
                maxID <- maxID + gene.num
                
                # all other attributes
                map <- data.frame(s=which(m.sel) + gene.num -1,
                      e=c((which(m.sel) + gene.num -1)[-1]-1,nrow(tmp.gff)),
                      g=queryHits(findOverlaps(g.grng,m.grng)),
                      o=oldmID)
                
                # mRNA IDs
                mids <- lapply(elementNROWS(split(map[,"g"],map[,"g"])),":",1)
                if(tmp.gff[1,]$strand == "+"){
                  mids <- lapply(mids,rev)
                }
                map$m=unlist(mids,use.names=FALSE)
                
                tmp.gff[unlist(mapply(":",map[,"s"],map[,"e"])),]$gffAttributes <- 
                  unlist(lapply(1:nrow(map),function(i){
                    ro <- map[i,]
                    m <- paste0(newIDs[ro$g],".",ro$m)
                    c(paste0("ID=",m,";Parent=",newIDs[ro$g]),
                      gsub(sub(sprintf("%s_",seqnames(tmp.gff)[1]),"",ro$o),
                           m,tmp.gff$gffAttributes[(ro$s+1):ro$e]))
                  }))
                
                if(fix){
                  # fix cds IDs - we may at least hope that
                  # PASA is consistent enough in its naming
                  # and positioning.
                  lIDs <- getGffAttribute(tmp.gff,"ID")[,1]
                  c.sel <- grep("cds",lIDs)
                  e.sel <- c.sel - 1
                  
                  # simple check
                  stopifnot(length(grep("exon",lIDs[e.sel])) == length(e.sel))
                  tmp.gff[c.sel]$gffAttributes <- sub("exon","cds",tmp.gff[e.sel]$gffAttributes)
                  
                  lIDs <- getGffAttribute(tmp.gff,"ID")[,1]
                  if(any(duplicated(lIDs))){
                   message("Fixing duplicated feature IDs (non mRNA or gene)")
                    dIDs <- lIDs[duplicated(lIDs)]
                    dIDs.sel <- lIDs %in% dIDs
                    # check
                    stopifnot(all(unique(tmp.gff[dIDs.sel]$type) != c("gene","mRNA")))
                    
                    # quick fix, just make them unique
                    tmp.gff$gffAttributes[dIDs.sel] <-
                      sprintf("ID=%s;%s",
                              make.unique(lIDs[dIDs.sel]),
                              sub(".*;Parent","Parent",tmp.gff$gffAttributes[dIDs.sel]))
                  }
                }
                
                # and extend
                split.gff <- c(tmp.gff,split.gff)
                rm(tmp.gff)
                
                # define the ranges to remove; get until the next gene
                se.rm <- c(se.rm,spanPos)
                               
                # create a new df
                prev <- prev[rep(1,gene.num),]
                prev$newID <- newIDs
                df.split <- rbind(df.split,prev)
              }
            }
            
            # update the df
            df.pos <- match(df[df$status=="split","oldID"],df.split$oldID)
            df[df$status=="split",] <- df.split[df.pos,]
            df.split <- df.split[-df.pos,]
            
            ## ------------
            ## GENE MERGE
            ## ------------
            message(sprintf("Considering gene merge 4/%s",totalSteps))

            # update the gene IDs for merged and novel
            newIDs.sel <- !grepl("^[A-Z][a-z]{4}[0-9]{6}g[0-9]{5}$",IDs[sel])
            nov.sel <- grepl("novel_gene",IDs[sel][newIDs.sel]) & df$status[newIDs.sel] != "split"
            
            # merged genes
            df$status[newIDs.sel][!nov.sel] <- "merged"
            
            # identify if there are multi merge
            gene.rm=NULL
            multi.merge <- lapply(strsplit(df$oldID[newIDs.sel][!nov.sel],"_"),"%in%",IDs[sel])
            
            # create the new ID
            df$newID[newIDs.sel][!nov.sel] <- sub("_.*","",df$oldID[newIDs.sel][!nov.sel])
            
            if(any(sapply(multi.merge,sum) > 0)){
              for (i in which(sapply(multi.merge,any))){
                
                ovl.id <- strsplit(df[newIDs.sel,][!nov.sel,][i,"oldID"],"_")[[1]][multi.merge[[i]]]
                #message(ovl.id)
                
                ## if the merged ID is different use this one (it's easier, less gff editing)
                df$newID[newIDs.sel][!nov.sel][i] <- ovl.id
                
                ## NOT IMPLEMENTED AS A LOOP; FIX IF IT HAPPENS
                stopifnot(length(ovl.id)==1)
                
                # get the range (we need the full range because of the previous sorting that
                # mixed the 2 genes, since their coordinates overlap)
                ovl.gene <- which(sel)[which(IDs[sel] == ovl.id)]
                kpt.gene <- which(sel)[which(IDs[sel] == df[newIDs.sel,][!nov.sel,][i,"oldID"])]
                start <- min(ovl.gene,
                             kpt.gene)
                end <- max(which(sel)[which(IDs[sel] == ovl.id)+1] - 1,
                           which(sel)[which(IDs[sel] == df[newIDs.sel,][!nov.sel,][i,"oldID"])+1] - 1)
                
                ## now edit the attribute
                gff.attrs <- gff[start:end,]$gffAttributes
                gff.attrs[grep("_[A-Z][a-z]{4}",
                               gff[start:end,]$gffAttributes)] <- gsub("[A-Z][a-z]{4}[0-9]+[^;]+g[0-9]+",
                                                                       ovl.id,
                                                                       gff.attrs[grep("_[A-Z][a-z]{4}",
                                                                                      gff[start:end,]$gffAttributes)])
                gff[start:end,]$gffAttributes <- gff.attrs
                
                ## update the gene coordinates to the max ones
                gff[kpt.gene,1:2] <- range(gff[c(ovl.gene,kpt.gene),1:2])
                
                ## check for duplicated IDs (escaping gene and CDS)
                ## there should no be any dups at that stage
                stopifnot(!any(duplicated(gsub(";.*","",gff.attrs[! gff[start:end,]$type  %in% c("gene","CDS")]))))
                
                ## update the status
                df$status[df$oldID == ovl.id] <- "merged.as.sv"
              
                ## record the gene
                gene.rm=c(gene.rm,which(sel)[which(IDs[sel] == ovl.id)])
              }
            }
            
            # create a df extension
            df.addID <- do.call(rbind,lapply(strsplit(df$oldID[newIDs.sel][!nov.sel],"_"),function(l){
              data.frame(oldID=l[2:length(l)],
                         newID=rep(l[1],length(l)-1),
                         status="merged",
                         scaffold=sub("g.*","",l[1]),
                         stringsAsFactors=FALSE)
            }))
            
            ## remove the genes edited previously (multi-merge)
            df.addID <- df.addID[ !df.addID$newID %in% gene.rm,]
            
            ## edit the gff to change all the IDs
            ## We edit everything on the same scaffold
            ## It is "a bit" overdoing it, but at least, it should
            ## avoid issues with the sorting
            m.pos <- do.call(rbind,mclapply(df[newIDs.sel,"scaffold"][!nov.sel],
                                       function(scf,gff){range(which(seqnames(gff) == scf))},
                                       gff,mc.cores=mc.cores))
            pos <- mapply(":",m.pos[,1],m.pos[,2])  
            gff[unlist(pos),]$gffAttributes <- gsub("_[A-Z][^;]*","",
                                                    gsub("_[A-Z][^;]*\\.[0-9]+","",
                                                         gff[unlist(pos),]$gffAttributes))
            
            ## ------------
            ## NOVEL GENES
            ## ------------
            message(sprintf("Creating novel gene IDs 5/%s",totalSteps))
            df$status[newIDs.sel][nov.sel] <- "novel"
            
            df[df$status=="novel","newID"] <- sprintf("%sg%05d",
                                                      df[df$status=="novel","scaffold"],
                                                      maxID + 1:sum(df$status=="novel"))
            
            ## DUPLICATES
            message(sprintf("Considering duplicate novel gene 6/%s",totalSteps))
            inx <- which(sel)[match(df$oldID[newIDs.sel][nov.sel],IDs[sel])]
            key <- paste(seqnames(gff[inx]),gff[,1][inx],gff[,2][inx],sep="_")
            dup.key <- duplicated(key)
            
            ## we only consider one duplicate
            stopifnot(all(table(key[dup.key]) == 1))
            
            if(sum(dup.key)>0){
              for(k in key[dup.key]){
                df.pos <- which(key == k)
                
                ## check the strand
                if(strand(gff[inx[df.pos[1]],]) == strand(gff[inx[df.pos[2]],])){
                  
                  ## same strand - proceed
                  ## to prevent issues further on, we flag the
                  ## first gene occurence as duplicate
                  df$status[newIDs.sel][nov.sel][df.pos[1]] <- "duplicate"
                  df$newID[newIDs.sel][nov.sel][df.pos[2]] <- df$newID[newIDs.sel][nov.sel][df.pos[1]]
                  df$newID[newIDs.sel][nov.sel][df.pos[1]] <- df$oldID[newIDs.sel][nov.sel][df.pos[2]]
                  
                  ## now update the geneID in the gff
                  m.sel <- which(Parents %in% df$oldID[newIDs.sel][nov.sel][df.pos])
                  p.sel <- which(Parents %in% IDs[m.sel])
                  spanPos <- sort(c(m.sel,p.sel))
                  
                  # all IDs/Parents
                  IDsParents <- as.data.frame(getGffAttribute(gff[spanPos],c("ID","Parent")),
                                              stringsAsFactors=FALSE)
                  
                  # update the mRNA Parent
                  IDsParents[spanPos %in% m.sel,"Parent"] <- df$newID[newIDs.sel][nov.sel][df.pos[2]]
                  
                  # update the feats
                  parts <- sub(paste(IDsParents[spanPos %in% m.sel,"ID"],collapse="|"),"",
                      IDsParents[spanPos %in% p.sel,"ID"])
                  
                  c.sel <- grep("cds",parts)
                  e.sel <- c.sel -1 
                  ## check further down in the code if that breaks.
                  stopifnot(length(grep("exon",parts[e.sel])) == length(e.sel))
                  parts[c.sel] <- sub("exon","cds",parts[e.sel])
                  
                  IDsParents[spanPos %in% p.sel,"ID"] <- 
                    sprintf("%s.%i%s",
                            df$newID[newIDs.sel][nov.sel][df.pos[2]],
                            rowSums(sapply(m.sel,"<",p.sel)),
                            parts)
                  
                  IDsParents[spanPos %in% p.sel,"Parent"] <- 
                    sprintf("%s.%i",
                            df$newID[newIDs.sel][nov.sel][df.pos[2]],
                            rowSums(sapply(m.sel,"<",p.sel)))
                  
                  # update the mRNA ID
                  IDsParents[spanPos %in% m.sel,"ID"] <- 
                    sprintf("%s.%i",
                            IDsParents[spanPos %in% m.sel,"Parent"],
                            1:length(m.sel))
                  
                  # and the attrs
                  gff[spanPos,]$gffAttributes <- 
                    sprintf("ID=%s;Parent=%s",IDsParents[,"ID"],IDsParents[,"Parent"])

                  ## and extend the gene.rm
                  gene.rm <- c(gene.rm,inx[df.pos[1]])
                  
                } else {
                  # different strand
                  df$status[newIDs.sel][nov.sel][df.pos[2]] <- "antisense"
                  
                  ## now update the geneID in the gff
                  m.sel <- which(Parents %in% df$oldID[newIDs.sel][nov.sel][df.pos[2]])
                  g.sel <- which(IDs %in% Parents[m.sel])
                  p.sel <- which(Parents %in% IDs[m.sel])
                  spanPos <- sort(c(m.sel,p.sel))
                  gff[g.sel]$gffAttributes <- sprintf("ID=%s;Name=%s",
                                                     df$newID[newIDs.sel][nov.sel][df.pos[2]],
                                                     df$newID[newIDs.sel][nov.sel][df.pos[2]])
                  
                  # all IDs/Parents
                  IDsParents <- as.data.frame(getGffAttribute(gff[spanPos],c("ID","Parent")),
                                              stringsAsFactors=FALSE)
                  
                  # update the mRNA Parent
                  IDsParents[spanPos %in% m.sel,"Parent"] <- df$newID[newIDs.sel][nov.sel][df.pos[2]]
                  
                  # update the feats
                  parts <- sub(paste(IDsParents[spanPos %in% m.sel,"ID"],collapse="|"),"",
                               IDsParents[spanPos %in% p.sel,"ID"])
                  
                  c.sel <- grep("cds",parts)
                  e.sel <- c.sel -1 
                  ## check further down in the code if that breaks.
                  stopifnot(length(grep("exon",parts[e.sel])) == length(e.sel))
                  parts[c.sel] <- sub("exon","cds",parts[e.sel])
                  
                  IDsParents[spanPos %in% p.sel,"ID"] <- 
                    sprintf("%s.%i%s",
                            df$newID[newIDs.sel][nov.sel][df.pos[2]],
                            rowSums(sapply(m.sel,"<",p.sel)),
                            parts)
                  
                  IDsParents[spanPos %in% p.sel,"Parent"] <- 
                    sprintf("%s.%i",
                            df$newID[newIDs.sel][nov.sel][df.pos[2]],
                            rowSums(sapply(m.sel,"<",p.sel)))
                  
                  # update the mRNA ID
                  IDsParents[spanPos %in% m.sel,"ID"] <- 
                    sprintf("%s.%i",
                            IDsParents[spanPos %in% m.sel,"Parent"],
                            1:length(m.sel))
                  
                  # and the attrs
                  gff[spanPos,]$gffAttributes <- 
                    sprintf("ID=%s;Parent=%s",IDsParents[,"ID"],IDsParents[,"Parent"])
                }
              }
              ## update the nov.sel selector to remove the duplicatess
              nov.sel[which(nov.sel)[df$status[newIDs.sel][nov.sel] %in% 
                                       c("antisense","duplicate")]] <- FALSE
            }
            
            ## NOVEL
            message(sprintf("Considering novel gene 7/%s",totalSteps))
            
            ## update the gff
            nov.genes <- df[newIDs.sel,][nov.sel,]
            posAttr <- do.call(rbind,mclapply(1:nrow(nov.genes),function(i,n){
            
              # get the entries spanned by that gene
              g.id <- which(IDs %in% n[i,"oldID"])
              m.id <- which(Parents %in% n[i,"oldID"])
              p.id <- which(Parents %in% IDs[Parents %in% n[i,"oldID"]])
              gffPos <- sort(c(g.id,m.id,p.id))
              
              # attrs
              tmp.attrs <- gff[gffPos]$gffAttributes
              
              # update the gene name
              tmp.attrs[which(gffPos %in% g.id)] <- 
                    sprintf("ID=%s;Name=%s",
                            n[i,"newID"],
                            n[i,"newID"])

              # update the mRNA
              tmp.attrs[which(gffPos %in% m.id)] <- 
                     sprintf("ID=%s.%i;Parent=%s",
                             n[i,"newID"],
                             1:length(m.id),
                             n[i,"newID"])
              
              # update the feats
              parts <- gsub("ID=cds.|ID=|;Parent=","",
                            gsub(paste0(IDs[m.id],collapse="|"),
                                 "",tmp.attrs[which(gffPos %in% p.id)]))
              part.sel <- which(parts=="")
              if(length(part.sel)){
                stopifnot(length(grep("exon",parts[part.sel-1]))==length(part.sel))
                parts[part.sel] <- sub("exon","cds",parts[part.sel-1])
              }
              tmp.attrs[which(gffPos %in% p.id)] <- 
                sprintf("ID=%s.%i%s;Parent=%s.%i",
                        n[i,"newID"],
                        rowSums(sapply(m.id,"<",p.id )),
                        parts,
                        n[i,"newID"],
                        rowSums(sapply(m.id,"<",p.id ))
                )
              
              # report
              data.frame(pos=gffPos,
                         atr=tmp.attrs,stringsAsFactors = FALSE)
            },nov.genes,mc.cores=mc.cores))
              
            # update the attrs
            gff[posAttr$pos]$gffAttributes <- posAttr$atr
            
            ## ------------
            ## FINALISE
            ## ------------
            message(sprintf("Finalizing 8/%s",totalSteps))
            
            # update the source
            if(!is.na(src)){
              levels(gff$source) <- c(levels(gff$source),src)
              gff$source[is.na(gff$source)] <- src
            }
            
            # update the gff
            ## create  list of all entry index to be removed
            ## and remove them all at once
            if(!is.null(se.rm)){
              gene.rm <- sort(unique(c(gene.rm,se.rm)))
            }
            if(!is.null(gene.rm)){
              gff <- gff[-gene.rm,]
            }
            
            if(!is.null(split.gff)){
              gff <- c(gff,split.gff)
            }
            
            # update the df
            if(!is.null(df.split)){
              df <- rbind(df,df.split)
            }
            df <- rbind(df,df.addID)
            
            # check duplicated IDs and fix if required
            if(fix){
              IDsParents <- getGffAttribute(gff,c("ID","Parent"))
              
              dup.sel <- IDsParents[,"ID"] %in% IDsParents[,"ID"][duplicated(IDsParents[,"ID"])]
              
              ## all CDS
              c.sel <- which(dup.sel & gff$type=="CDS")
              e.sel <- c.sel -1 
              
              if(any(gff$type[e.sel]!="exon")){
                warning(sprintf("The mRNA %s has improperly formatted feature IDs\n",
                               unique(getGffAttribute(gff[e.sel][gff[e.sel]$type != "exon"],"Parent"))))
              }
              
              gff[c.sel]$gffAttributes <- sub("exon","cds",gff[e.sel]$gffAttributes)
            
              ## the rest
              dup.sel <- which(dup.sel & gff$type!="CDS")
              attr <- gff[dup.sel]$gffAttributes
              
              ## no gene
              stopifnot(table(gff[dup.sel]$type)["gene"] == 0)
              warning(sprintf("%i more feature IDs are improperly formatted",length(dup.sel)))
              
              # all IDs
              IDsParents[dup.sel,"ID"] <- make.unique(IDsParents[dup.sel,"ID"])
              
              ## mRNA as Parents
              m.sel <- which(gff[dup.sel]$type == "mRNA")
              IDsParents[dup.sel,"Parent"][-m.sel] <- rep(IDsParents[dup.sel,"ID"][m.sel],
                                                          diff(c(m.sel,length(dup.sel)+1))-1)
              gff[dup.sel]$gffAttributes <- sprintf("ID=%s;Parent=%s",
                                                    IDsParents[dup.sel,"ID"],
                                                    IDsParents[dup.sel,"Parent"])
            }
            
            # final gff sort
            #gff <- gff[order(seqnames(gff),
            #                    gff[,1],
            #                    gff$type
            #),]
            
            # write it to the output
            writeGff3(gff,file = output)
            
            # return the df ordered by newIDs
            return(df[order(df$newID),])
            
          })

#' ## validateProcessedGff3
#' This function takes three arguments:
#' 
#' 1. gff3: the input file - mandatory
#' 
#' 2. types: the gff type to check for. Defaults to c("gene","mRNA","exon","five_prime_UTR","three_prime_UTR")
#' 
#' 3. hierarchy: a matrix defining the relationship Parent <-> ID of the above types
#' 
#' This function invisibly return TRUE if the validation succeed. It returns
#' FALSE otherwise and warns about the encountered issues
setMethod(f = "validateProcessedGff3", 
          signature = "character",
          definition = function(
            gff3=character(0),
            types=c("gene","mRNA","exon","five_prime_UTR","three_prime_UTR","CDS"),
            hierarchy=matrix(c("gene",rep("mRNA",5),"exon",
                               "five_prime_UTR","three_prime_UTR","CDS"),ncol=2)
          ){

            # some checks
            stopifnot(length(gff3)>0)
            stopifnot(file.exists(gff3))
            stopifnot(length(types) > 0)
            
            # then proceed
            message("Reading the gff3 file 1/4")
            gff <- readGff3(gff3,quiet=TRUE)
            
            # one more check
            stopifnot(any(levels(gff$type) %in% types))
            
            # check that IDs are unique
            res <- sapply(types,function(typ){
              message(paste("Checking IDs of type:",typ,"2/4"))
              dup <- duplicated(getGffAttribute(gff[gff$type==typ,],"ID"))
              if(any(dup)){
                warning(sprintf("Attribute of type %s with ID %s is duplicated\n",
                                typ,getGffAttribute(gff[gff$type==typ,],"ID")[dup]))
              }
              sum(dup)
            })
            
            # check that IDs do not overlap groups
            message("Checking all IDs 3/4")
            dups <- duplicated(getGffAttribute(gff[gff$type %in% types,],"ID"))
            if(sum(dups > 0)){
              warning(sprintf("The ID %s is duplicated in an attribute of type %s\n",
                getGffAttribute(gff[gff$type %in% types,],"ID")[dups],
                gff[gff$type %in% types,]$type[dups]))
            }
            res <- c(res,sum(dups))
            
            # check that Parents all exists
            message("Checking Parent - ID relationships 4/4")
            res <- c(res,apply(hierarchy,1,function(ro){
              IDs <- getGffAttribute(gff[gff$type == ro[[1]],],"ID")
              # message(sprintf("We have %s IDs",length(IDs)))
              Parents <- getGffAttribute(gff[gff$type == ro[[2]],],"Parent")
              # message(sprintf("We have %s Parents",length(Parents)))
              fail <- Parents[!Parents %in% IDs]
              # message(sprintf("And %s orphans",length(fail)))
              if(length(fail) > 0){
                warning(sprintf("The %s ID %s referenced as a Parent does not exist\n",
                                ro[[1]],fail))
              }
              return(length(fail))
            }))            
            
            # done
            invisible(ifelse(sum(res)==0,TRUE,FALSE))
          })

#' # Caveat
#' There are a number of caveats, since we make a large number of assumptions
#' when processing the IDs.
#' Namely:
#' 
#' 1. we expect IDs of the form [A-Z][a-z]{4}[0-9]{6}g[0-9]{5}
#' 
#' 2. we expect the gene IDs to be sorted by start position and hence based
#' the ID renumbering on that.
#' 
#' It is essential that after running the **checkAndUpdatePasaGff3IDs**, you
#' validate the obtained gff using the **validateProcessedGff3** function.
#' 
