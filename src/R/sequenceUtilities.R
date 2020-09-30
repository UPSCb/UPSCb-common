#' ---
#' title: "Sequence utilities"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions that help extract
#' various sequences from a genome using a gff3 file
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/sequenceUtilities.R")
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
stopifnot(suppressPackageStartupMessages(require(Biostrings)))
stopifnot(suppressPackageStartupMessages(require(genomeIntervals)))

#' 
#' # Generics
setGeneric(name="extract",
           def=function(gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
                        genome=DNAStringSet(),
                        feature=c("CDS","mRNA","protein")
           ){
             standardGeneric("extract")
           })

setGeneric(name="extractFromGenome",
           def=function(gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
                        genome=DNAStringSet()
                        ){
             standardGeneric("extractFromGenome")
           })

setGeneric(name="extractCdsFromGenome",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             genome=DNAStringSet()
           ){
             standardGeneric("extractCdsFromGenome")
           })

setGeneric(name="extractMrnaFromGenome",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             genome=DNAStringSet()
           ){
             standardGeneric("extractMrnaFromGenome")
           })

setGeneric(name="extractProteinFromGenome",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             genome=DNAStringSet()
           ){
             standardGeneric("extractProteinFromGenome")
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
setMethod(f = "extract", 
          signature = c("Genome_intervals","DNAStringSet","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            genome=DNAStringSet(),
            feature=c("CDS","mRNA","protein")
          ){
            
            ## validation
            if(nrow(gff3) == 1 & gff3[1,1] == 0 & gff3[1,2] == 0){
              stop("You need to provide a non empty Genome_intervals_stranded 'gff3' object")
            }            
            stopifnot(length(genome)>0)            
            feature <- match.arg(feature)
            
            ## dispatch
            switch(feature,
                   "CDS" = extractCdsFromGenome(gff3,genome),
                   "mRNA" = extractMrnaFromGenome(gff3,genome),
                   "protein" = extractProteinFromGenome(gff3,genome))
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
            
            # sanity check
            stopifnot(!all(duplicated(getGffAttribute(gff3,"ID")[,1])))
            
            # Get all sub sequences
            sseq <- subseq(genome[seqnames(gff3)],gff3[,1],gff3[,2])
            
            # Get the parent
            prt <- getGffAttribute(gff3,"Parent")
            
            # check for the absence of parent (e.g. type == gene)
            if(all(is.na(prt))){
              prt <- getGffAttribute(gff3,"ID")
            }

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