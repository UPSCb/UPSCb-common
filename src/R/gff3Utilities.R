#' ---
#' title: "Gff3 utilities"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # What is this about?
#' 
#' This file is only a source file containing functions that help subtract
#' a gff3 file by a list of IDs
#' 
#' To use it in your script do:
#' ```{r, eval=FALSE}
#' source("~/Git/UPSCb/src/R/gff3Utilities.R")
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

#' 
#' # Generics
setGeneric(name="extract",
           def=function(gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
                        IDs=character(0),
                        feature=c("gene","mRNA","exon")
           ){
             standardGeneric("extract")
           })

setGeneric(name="extractFromGff3UsingGeneIDs",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             IDs=character(0)
           ){
             standardGeneric("extractFromGff3UsingGeneIDs")
           })

setGeneric(name="extractFromGff3UsingMrnaIDs",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             IDs=character(0)
           ){
             standardGeneric("extractFromGff3UsingMrnaIDs")
           })

setGeneric(name="extractFromGff3UsingExonIDs",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             IDs=character(0)
           ){
             standardGeneric("extractFromGff3UsingExonIDs")
           })

setGeneric(name="getIDPositions",
           def=function(
             gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
             IDs=character(0),
             feature=c("gene","mRNA","exon")
             ){
             standardGeneric("getIDPositions")
           })

#' # Implementation
#' 
#' ## Arguments
#' The 'extract' function share the following arguments:
#' 
#' 1. gff3: a gff3 read as a Genome_intervals_stranded object using the genomeIntervals package readGff3 function.
#' 
#' 2. IDs: a list of IDs (as a character vector) for which the gff3 structure needs to be extracted.
#'
#' The extract function (the most generic) has a third argument to define the type of feature to retrieve
#' 
#' 3. feature: one of 'gene', 'mRNA' or 'exon'. 
#' This represent the hierarchy of a gff3 file that will be traversed. Hence, any
#' feature whose ID is in the IDs list or a child thereof 
#' (i.e. their Parent attributes is used to retrieve their parent at the first and or second level) 
#' or inheriting from it (they are parent of other features are the first or second level) will
#' be retrieved irrespective of the gff3 feature type' e.g. CDS features will also be retrieved. 
#'
#' The 'getIDPositions' function takes the same arguments as extract, but should be considered internal
#' and not called directly.
#'
#' ## Functions
#' ### extract
#'  
#' This function is a wrapper that based on the desired feature delegates
#' to the appropriate "worker".
setMethod(f = "extract", 
          signature = c("Genome_intervals","character","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            IDs=character(0),
            feature=c("gene","mRNA","exon")
          ){

            ## dispatch
            switch(match.arg(feature),
                   "gene" = extractFromGff3UsingGeneIDs(gff3,IDs),
                   "mRNA" = extractFromGff3UsingMrnaIDs(gff3,IDs),
                   "exon" = extractFromGff3UsingExonIDs(gff3,IDs))
          })

#' ### extractFromGff3UsingGeneIDs
#' This function is one of the "worker" that extract the gene feature and
#' descendants from a gff3 file based on the provided gene IDs
#' 
setMethod(f = "extractFromGff3UsingGeneIDs", 
          signature = c("Genome_intervals","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            IDs=character(0)
          ){
            gff3[getIDPositions(gff3,IDs,"gene"),]
          })

#' ### extractMrnaFromGenome
#' This function is one of the "worker" that extract the gene feature and
#' descendants from a gff3 file based on the provided mRNA IDs
#' 
setMethod(f = "extractFromGff3UsingMrnaIDs", 
          signature = c("Genome_intervals","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            IDs=character(0)
          ){
            gff3[getIDPositions(gff3,IDs,"mRNA"),]
          })

#' ### extractProteinFromGenome
#' This function is one of the "worker" that extract the gene feature and
#' descendants from a gff3 file based on the provided exon IDs
#' 
setMethod(f = "extractFromGff3UsingExonIDs", 
          signature = c("Genome_intervals","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            IDs=character(0)
          ){
            gff3[getIDPositions(gff3,IDs,"exon"),]
          })

#' ### getIDPositions
#' This function returns an integer vector of the position of the
#' feature of interest in the gff3 file. This ideally should be an internal
#' file as it contains the ascendants/descendants logic
setMethod(f = "getIDPositions",
          signature = c("Genome_intervals","character","character"),
          definition = function(
            gff3=GenomeIntervals(chromosome="empty",start=0,end=0),
            IDs=character(0),
            feature=c("gene","mRNA","exon")
            ){
            
            ## validation
            if(nrow(gff3) == 1 & gff3[1,1] == 0 & gff3[1,2] == 0){
              stop("You need to provide a non empty Genome_intervals_stranded 'gff3' object")
            }
            stopifnot(length(IDs)>0)
            feature <- match.arg(feature)
            
            ## get the ID and Parent
            ID <- getGffAttribute(gff3,"ID")
            Parent <- getGffAttribute(gff3,"Parent")
            
            ## get the positions
            return(.getPosition(feature,ID,Parent,IDs))
          })

#' # An internal function
".getPosition" <- function(feature,ID,Parent,IDs){
  sort(switch(feature,
              "gene" = {c(which(ID %in% IDs),
                          which(Parent %in% IDs),
                          which(Parent %in% ID[Parent %in% IDs]))},
              "mRNA" = {c(which(ID %in% Parent[ID %in% IDs]),
                          which(ID %in% IDs),
                          which(Parent %in% IDs))},
              "exon" = {c(which(ID %in% Parent[ID %in% Parent[ID %in% IDs]]),
                          which(ID %in% Parent[ID %in% IDs]),
                          which(ID %in% IDs))}))
}

#' # SessionInfo
#' ```{r for the doc}
#' sessionInfo()
#' ```
