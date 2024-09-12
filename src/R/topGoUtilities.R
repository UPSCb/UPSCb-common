#' ---
#' title: "topGO at UPSC"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Disclaimer
#' 
#' **1. This function has hardcoded paths valid only for the UPSC.**
#' 
#' **2. This function has no input check whatsoever at present, use with care**
#'
#' # TODO
#' 
#' 1. Fix 2. above.
#' 2. Consider integrating more methods and their comparisons through the GenTable function
#'
#' # Libraries
suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(topGO)
}))

#' 
#' # Gene to GO mappings
#' 
#' This function retrieves existing mappings. By default it retrieved mappings listed in a file, unless `refresh=TRUE`
#' 
listMappings <- function(refresh=FALSE){
  if(refresh){
    list.files("/mnt/picea/storage/reference",pattern="gene_to_go.tsv",full.names=TRUE,recursive=TRUE)    
  } else {
    scan("/mnt/picea/storage/reference/topGO/mappingList_20231205.txt",what="character")
  }
}

#' # Prepare the selected annotation
#' 
#' This function takes one of the mapping provided above and converts them into the object required by topGO
prepAnnot <- function(mapping){
  annot <- read_delim(mapping,col_names=FALSE,show_col_types=FALSE)
  geneID2GO <- lapply(lapply(unlist(annot[,2],use.names=FALSE),strsplit,"\\|"),unlist)
  names(geneID2GO) <- unlist(annot[,1],use.names=FALSE)
  return(geneID2GO)
}

#' # Perform the enrichment
#' 
#' This function takes as mandatory input:
#' 1. gene set - list of gene IDs of "interest"
#' 2. background - list of gene IDs in the population
#' 3. annotation - the gene to GO mapping created by the prepAnnot function.
#' 
#' Optional arguments:
#' 1. All three ontologies are searched by default
#' 2. The default algorithm is parent-child
#' 3. The corresponding test is fisher
#' 4. The default filtering is FDR based (Benjamini-Hochberg correction), at a 1% cutoff. Both method and significance cutoff (alpha) can be controlled. Check `p.adjust.methods` for possible choices of correction including "none".
#' 
#' Note: Check the topGO vignette for alternative algorithms / statistics and their possible combination (see Table 1 in the Introduction section)
#' 
#' Value: it returns a list of length "ontology" of tibbles with six columns. Names are self-explanatory.
#'  
topGO <- function(set,background,annotation,
                  ontology=c("BP","CC","MF"),
                  algorithm="parentchild",
                  statistic="fisher",
                  p.adjust=sort(p.adjust.methods),
                  alpha=0.01,
                  getgenes=FALSE){
  
  p.adjust <- match.arg(p.adjust)
  
  # create the allGenes
  allGenes <- factor(as.integer(background %in% set))
  names(allGenes) <- background
  
  # iterate over the ontologies
  lst <- lapply(ontology, function(o,g,a){
    GOdata <- new("topGOdata", 
                  ontology = o, 
                  allGenes = g,
                  annot = annFUN.gene2GO, 
                  gene2GO = a)
    allGO <- genesInTerm(GOdata)
    
    results <- runTest(GOdata,algorithm=algorithm,statistic=statistic)
    
    n <- ifelse(p.adjust=="none",
                sum(score(results) <= alpha),
                sum(p.adjust(score(results),method=p.adjust) <= alpha))

    if(n==0){
      return(NULL)
    } else{
      resultTable <- as_tibble(GenTable(GOdata,
                                        results,
                                        topNodes=n)) %>% 
                       rename_with(function(sel){"FDR"},.cols=last_col())
      if(getgenes){
        resultTable <- resultTable %>%
          mutate(allgenes = map(GO.ID,function(x){allGO[[x]]})) %>%
          mutate(siggenes = map(allgenes,function(x){unlist(x)[unlist(x) %in% set]}))
      }
      resultTable
    }
  },allGenes,annotation)
  names(lst) <- ontology
  return(lst)
}

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```
