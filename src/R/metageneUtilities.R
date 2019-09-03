#' ---
#' title: "Metagene utilities"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions that help 
#' (i) extract metagene bins from a gff3 file (createBins)
#' (ii) summarise coverage/expression values per metagene (computeBins)
#' (iii) plot the metagenes (plotMetagenes)
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/metageneUtilities.R")
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
stopifnot(suppressPackageStartupMessages(require(GenomicRanges)))

#' 
#' # Generics
setGeneric(name="createBins",
           def=function(gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
                        feature=c("exon"),
                        flank=1000L
           ){
             standardGeneric("createBins")
           })

setGeneric(name="computeBins",
           def=function(bins=list(),
                        genome=DNAStringSet()
           ){
             standardGeneric("extractFromGenome")
           })

setGeneric(name="getProteinState",
           def=function(
             prot=AAStringSet()
           ){
             standardGeneric("getProteinState")
           })

#' # Implementation
#' 
#' ## Arguments
#' The 'extract' function share the following arguments:
#' 
#' 1. gff3: a gff3 read as a Genome_intervals_stranded object using the genomeIntervals package readGff3 function.
#' 
#' 2. genome: a genome read as a DNAStringSet object using the Biostrings package readDNAStrinSet function.
#'
#' The extract function (the most generic) has a third argument to define the type of feature to retrieve
#' 
#' 3. feature: one of 'CDS', 'mRNA' or 'protein'
#'
#' The 'getProteinState' function takes an _AAStringSet_ as argument.
#'
#' ## Functions
#' ### extract
#'  
#' This function is a wrapper that based on the desired feature delegates
#' to the appropriate "worker".
setMethod(f = "createBins", 
          signature = c("Genome_intervals"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            feature=c("exon"),
            gene.nchunk=1000L,
            flank.bp=1000L,
            flank.nchunk=100L,
            mc.cores=1L
          ){
            
            ## validation
            if(nrow(gff3) == 1 & gff3[1,1] == 0 & gff3[1,2] == 0){
              stop("You need to provide a non empty Genome_intervals_stranded 'gff3' object")
            }
            
            ##
            mRNA.id.parent <- getGffAttribute(gff3[gff3$type=="mRNA",],c("ID","Parent"))
            
            ex.gff3 <- gff3[gff3$type==feature,]
            ex.grng <- as(ex.gff3,"GRanges")
            ex.glist <- split(ex.grng,mRNA.id.parent$Parent[match(as.character(ex.grng$Parent),mRNA.id.parent$ID)])
            r.ex.glist <- GenomicRanges::reduce(ex.glist)
            
            sums <- width(r.ex.glist)
            chunksize <- (sum(sums)-1)/(gene.nchunk)
            
            # TODO we need to do this by exon... i.e. assess the proportion of the 
            # exon and use that to define the number of windows
            
            GRangesList(mclapply(1:length(chunksize),
                     function(i,rgl,chk){
                       message(i)
                       ans_width <- diff(c(0,as.integer(
                         cumsum(rep.int(chk[i],gene.nchunk)))))
                       GRanges(seqnames=seqnames(rgl[[i]]),
                               ranges=IRanges(end=start(rgl[[i]])+cumsum(ans_width),
                                              width=ans_width+1),
                               strand=strand(r.ex.glist[[i]]))
            },r.ex.glist,chunksize,mc.cores=mc.cores))
            
          })

#' ### extractCdsFromGenome
#' This function is one of the "worker" that select the CDS
#' coordinates from the gff3 file and retrieve the corresponding sequence 
#' 
setMethod(f = "extractCdsFromGenome", 
          signature = c("Genome_intervals","DNAStringSet"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            genome=DNAStringSet()
          ){
            extractFromGenome(gff3[gff3$type=="CDS",],genome)
          })

#' ### extractMrnaFromGenome
#' This function is one of the "worker" that select the mRNA
#' coordinates from the gff3 file and retrieve the corresponding sequence 
#' 
setMethod(f = "extractMrnaFromGenome", 
          signature = c("Genome_intervals","DNAStringSet"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            genome=DNAStringSet()
          ){
            extractFromGenome(gff3[gff3$type=="exon",],genome)
          })

#' ### extractProteinFromGenome
#' This function is one of the "worker" that select the CDS
#' coordinates from the gff3 file and retrieve the corresponding 
#' protein sequence 
#' 
setMethod(f = "extractProteinFromGenome", 
          signature = c("Genome_intervals_stranded","DNAStringSet"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            genome=DNAStringSet()
          ){
            translate(extractFromGenome(
              gff3[gff3$type=="CDS",],
              genome), if.fuzzy.codon="solve")
          })

#' ### extractFromGenome
#' This function is the one having the common logic. It should be made
#' internal if that is ever packaged. Hence, **do no use it directly**
#' _unless you know what you are doing_.
#' 
setMethod(f = "extractFromGenome", 
          signature = c("Genome_intervals","DNAStringSet"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            genome=DNAStringSet()
          ){
            # order the gff by chrom, start, stop
            gff3 <- gff3[order(seqnames(gff3),gff3[,1]),]
            
            # Get all sub sequences
            sseq <- subseq(genome[seqnames(gff3)],gff3[,1],gff3[,2])
            
            # Get the parent
            prt <- getGffAttribute(gff3,"Parent")
            
            # Sort them
            DSSL <- split(sseq,prt)
            
            #' and collpase them into one sequence
            seq <- DNAStringSet(mclapply(1:length(DSSL),
                                         function(i,d){unlist(d[[i]])},
                                         DSSL,mc.cores=16))
            #' name all the sequences
            names(seq) <- names(DSSL)
            
            #' Reverse complement the - strand one
            strands <- gff3[match(names(seq),prt),]$strand
            seq[strands=="-"] <- 
              DNAStringSet(mclapply(which(strands=="-"),
                                    function(i,d){
                                      reverseComplement(d[[i]])
                                    },seq,mc.cores=16))
            
            #' return the sequences
            return(seq)
          })

#' ### getProteinState
#' This function returns a factor defining the protein state; i.e.
#' its coding potential: one of 'full-length','5p-partial','fragment',
#' '3p-partial' or 'non-coding'
setMethod(f = "getProteinState",
          signature = "AAStringSet",
          definition = function(
            prot=AAStringSet()){
            stp <- alphabetFrequency(prot)[,"*"]
            fp <- narrow(prot,start=1,width=1) == "M"
            tp <- subseq(prot,start=width(prot),width=1) == "*"
            
            return(factor(ifelse(fp,
                                 ifelse(stp==0,"5p-partial",
                                        ifelse(stp==1 & tp,"full-length",
                                               "non-coding")),
                                 ifelse(stp==0,"fragment",
                                        ifelse(stp==1 & tp,"3p-partial",
                                               "non-coding"))),
                          levels=c("full-length","5p-partial","fragment",
                                   "3p-partial","non-coding")))
          })