
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(purrr)
  library(topGO)
})

# Prepare annotation object
prepAnnot <- function(mapping){
  annot <- read_delim(mapping,col_names=FALSE,show_col_types=FALSE)
  geneID2GO <- lapply(lapply(unlist(annot[,2],use.names=FALSE),strsplit,"\\|"),unlist)
  names(geneID2GO) <- unlist(annot[,1],use.names=FALSE)
  return(geneID2GO)
}

# the following function gets a vector of genes with baseMean >= the baseMean of the gene with lower expression with non NA padj
get_background <- function(results_file) {
  suppressMessages(
    data <- read_csv(results_file, show_col_types = FALSE)
  )
  data <- data %>% dplyr::rename(gene = ...1)
  thres <- min(data[!is.na(data$padj), ]$baseMean)
  genes <- data %>% 
    filter(baseMean >= thres) %>% 
    dplyr::select(gene) %>% 
    pull()
  return(genes)
}


#' The following function is like the TopGO function found in the TopGO template, but
#' it returns a single dataframe with all enriched GO terms, rather than a list of three
#' dataframes with enriched biological processes, molecular functions and cellular components.

topGO_combined <- function(set,background,annotation,
                           ontology=c("BP","CC","MF"),
                           algorithm="parentchild",
                           statistic="fisher",
                           p.adjust=sort(p.adjust.methods),
                           alpha=0.05,
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
    results <- runTest(GOdata,algorithm=algorithm,statistic=statistic)
    
    n <- ifelse(p.adjust=="none",
                sum(score(results) <= alpha),
                sum(p.adjust(score(results),method=p.adjust) <= alpha))
    
    if(n==0){
      return(NULL)
    } else{
      #I added "numChar=1000 to not trim the GO
      resultTable <- as_tibble(GenTable(GOdata,
                                        results,
                                        numChar=1000,
                                        topNodes=n)) %>% 
        rename_with(function(sel){"FDR"},.cols=last_col())
      if(getgenes){
        allGO <- genesInTerm(GOdata)
        resultTable <- resultTable %>%
          mutate(allgenes = map(GO.ID,function(x){allGO[[x]]})) %>%
          mutate(siggenes = map(allgenes,function(x){unlist(x)[unlist(x) %in% set]}))  %>%
          mutate(allgenes = map(allgenes,function(x){paste(x,collapse = "|")}) %>% unlist(use.names = F),
                 siggenes = map(siggenes,function(x){paste(x,collapse = "|")}) %>% unlist(use.names = F))
      }
      resultTable
    }
  },allGenes,annotation)
  names(lst) <- ontology
  #This last bind_rows was added by me and is not part of the original function
  return(bind_rows(lst, .id = "GO category"))
}