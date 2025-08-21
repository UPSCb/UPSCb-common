## Doc (ideally roxygen syntax, us #' and keywords #' title)


## diff library vs require

suppressPackageStartupMessages({
  stopifnot({
    require(here)
    require(purrr)
    require(tidyverse)
    require(topGO)  
  })
})

## S4 functions
## 1. a generic
## 2. one or more method

# Prepare annotation object
## use key=value for the arguments 
prepAnnot <- function(mapping=character(0L)){
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
  
  ## require all packages
  data <- data %>% dplyr::rename(gene = ...1)
  thres <- min(data[!is.na(data$padj), ]$baseMean)
  data %<>% 
    filter(baseMean >= thres) %>% 
    dplyr::select(gene) %>% 
    pull()
  return(data)
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
  
  ## more validation of arguments
  
  # create the allGenes
  allGenes <- factor(as.integer(background %in% set))
  names(allGenes) <- background
  
  # iterate over the ontologies
  ## longer vaariable names (maybe...)
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
        rename_with(~"FDR",.cols=last_col())
      
      ## Always spell TRUE and FALSE
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
  ## NO space in column / variable names
  return(bind_rows(lst, .id = "GO_category"))
}