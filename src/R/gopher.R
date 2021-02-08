#' R wrapper for the gofer2 REST API
#'
#' @param genes 
#' @param background 
#' @param task 
#' @param alpha Maximun pAdj value permited
#' @param host 
#' @param port 
#' @param url 
#' @param endpoint 
#'
#' @return
#'
#' @examples
#' # Enrichment
#' enr <- gopher(c("K00001","K00002"),url='ko', task=list('ko_pathway'))
#' $ko_pathway
#' A tibble: 16 x 9
#' def                                             id           m    mt     n name        nt     padj      pval
#' <chr>                                           <chr>    <int> <int> <int> <chr>    <int>    <dbl>     <dbl>
#' 1 https://www.kegg.jp/dbget-bin/www_bget?map00010 map00010 22720   103     2 map00010     2 0.000326 0.0000204
#' 2 https://www.kegg.jp/dbget-bin/www_bget?map01220 map01220 22720   215     2 map01220     2 0.000713 0.0000891
#' 3 https://www.kegg.jp/dbget-bin/www_bget?map00626 map00626 22720    29     2 map00626     1 0.00516  0.00255  
#' 
#' # Term to gene lookup
#' Ex. y <- gopher(c("GO:0034641", "GO:0071705"), url="pabies", task=c("go"), endpoint = "term-to-gene")
#' y$`GO:0071705`
#' [1] "MA_10427660g0010" "MA_10436496g0010" "MA_552029g0010"   "MA_59139g0010"    "MA_847283g0010"  
#' y$`GO:0034641`
#' [1] "MA_10053307g0010" "MA_10263024g0010" "MA_10426126g0010" "MA_15844g0010"    "MA_169451g0010"   "MA_3184966g0010"  "MA_40229g0010"   
gopher <- function(genes=character(0),
                   background = NULL,
                   task = list("go",
                            "pfam",
                            "kegg",
                            "mapman"),
                   alpha=0.05,
                   host="https://franklin.upsc.se",
                   port=5432,
                   url="json",
                   endpoint=c("enrichment","gene-to-term","get-sets","term-to-gene"),
                   override=FALSE) {

  suppressPackageStartupMessages({
    require(tidyr)
    require(tibble)
    require(dplyr)
    require(jsonlite)
  })
  
  # arguments
  endpoint <- match.arg(endpoint)
  
  if(!endpoint %in% c("term-to-gene","gene-to-term")&& !is.list(task)){
    task <- as.list(task)
  }
  
  #term-to-gene behavior
  if(endpoint %in% c("term-to-gene","gene-to-term")) {
    if(length(task)>1) {
      stop("Term to gene, only accepts one type of enrichment per query.")
    }
    
    enrTerms <- genes
    # gofer2 doesn't admit 1 term only, a quick fix is to duplicate the same term,
    # it doesn't change the result.
    if (length(enrTerms) == 1) {
      enrTerms <- rep(enrTerms, 2)
    }
    # nullify background
    background = NULL
  }
  

  # sanity check
  # any transcript IDs?
  # CAVEAT: very UPSC specific (almost)
  # replace by a function that gets the complete pop from gopher and use that to test
  if(endpoint != "term-to-gene" && any(grepl("\\.\\d+$",genes))){
    stop("Your gene list contains transcript IDs (ending in .[0-9]+); trim them.")
  }
  
  
  if(!is.null(background)){
    if(length(background)==0){
      stop("Your background can not be empty. To provide no background (use the whole population), set it to NULL.")
    }
    # any transcript IDs?
    # CAVEAT: very UPSC specific (almost)
    # replace by a function that gets the complete pop from gopher and use that to test
    if(any(grepl("\\.\\d+$",background))){
      stop("Your gene list contains transcript IDs (ending in .[0-9]+); trim them.")
    }
    # all genes in background?
    if(!all(genes %in% background)){
      stop("Not all the genes you provide as a query are in your background")
    }
  }

  body <- switch (endpoint,
    "term-to-gene" = list(target = list(terms = enrTerms, name = task)),
    "gene-to-term" = list(genes = enrTerms, target = task),
    list(
      genes = genes,
      background = background,
      target = task,
      alpha = alpha
    ))
  
  # request
  # catch conection errors
  request <- tryCatch({
    httr::POST(paste0(host,":",port,"/",url,"/",endpoint), body=body, encode = "json")
    },
    error = function(e) {
      message(e)
      return(NULL)
    }
  )
  # if the connection fails return a NULL object
  if(is.null(request)){
    return(NULL)
  }
  
  # process
  # 
  parsed <- tryCatch({
    jsonlite::fromJSON(httr::content(request, as = "text", encoding = "UTF-8"))
    }, 
    error = function(e) {
      message(e)
      message("Response is not in a correct format.")
      message(request)
      return(NULL)
      }
    )
  
  # if parsing fails return a NULL object
  if(is.null(parsed))
    return(NULL)
  
  # term-to-gene parsing is different from enrichment
  # a list of terms will be returned, each list element will contain the name
  # of the genes related to that term.

  if(endpoint=="term-to-gene"){
    parsed <- parsed[[task]]
    termNames <- parsed$term
    parsed <- parsed[,"ids"]
    names(parsed) <- termNames
    #print(parsed)
    return(parsed)
  }
  
  if(endpoint=="gene-to-term"){
    return(tibble(ID=rep(parsed[[task]]$id,
                  sapply(parsed[[task]]$terms,nrow)),
           GOID=unlist(lapply(parsed[[task]]$terms,"[[","id"))))
  }
  np <- names(parsed)  
  
  if(!is.null(parsed$err)){
    warning(parsed$err)
    if(override==TRUE){
      np <- np[!grepl("err",np)]  
    } else {
      return(NULL)
    }
  }
  
  # return
  parsed <- lapply(np, function(n){
    f <- parsed[[n]]
  
    if(length(f) == 0 | is.null(f))
    {
      message("No enrichments found in task: ", n)
      return(NULL)
    }
    return(as_tibble(f) %>% arrange(padj))
  })
  names(parsed) <- np
  parsed <- parsed[unlist(task)]
  return(parsed)
}


#' function to collect and retrieve the most likely genes
#' that created an enrichment
#'
#' Note that because of the enrichment process, which can report GO terms enriched while they are
#' not associated to any of the genes, this function reverse engineer the most likely origin. It is
#' however possible that this process fails. The results contain an additional flag to report whether
#' the gene-term link was direct or the result of a walk in the GO graph. The reported association 
#' term - gene is then marked as ancestor to represent the fact that the provided GO term is an ancestor
#' of the term to which the gene is associated
#'
#' @param genes - a list of gene IDs (typically genes from a differential expression analysis)
#' @param terms - a list of GO terms (typically those resulting from the enrichment of the `genes`)
#' @param url - the gofer url to use, default to potra2
#' @param mc.cores - number of CPUs to use, default to 1
#'
#' @return a tibble with 3 columns: ID (the gene ID), GOID (the corresponding enriched term),
#' and link: the type of link, one of `direct`, `ancestor` or `missing`. `direct` means that the GOID was directly 
#' associated to a gene in genes, while `ancestor` signifies that the GO term reported is an 
#' ancestor of a term associated with a gene in genes. Finally, a third type of link can be reported:
#' `missing`, which signifies that either the term is obsolete (terms used an older version of GO) or
#' novel (the GO.db is outdated).
#'
#' @examples
#' genes <- c("Potra2n10c21792","Potra2n10c21793")
#' terms <- c("GO:0005886","GO:0016020","GO:0110165","GO:0051193")
#' enrichedTermToGenes(genes, terms)

enrichedTermToGenes <- function(genes,terms,url="potra2",mc.cores=1L){
  suppressPackageStartupMessages({
    require(GO.db)
    require(parallel)
    require(tidyverse)
  })
  
  # make sure that the terms are unique
  terms <- unique(terms)
  
  # get the search space, i.e. the GO terms associated with the genes
  associatedIDs <- gopher(genes,url=url,task=c("go"), endpoint = "gene-to-term")
  
  # find the terms either obsolete or missing in GO.db (the latter if GO.db is out-of-date)
  missing <- ! terms %in% keys(GO.db)
  
  # extract the direct annotation
  associatedIDs %>% filter(GOID %in% terms) %>% mutate(link="direct") %>%
    # report the missing
    bind_rows(tibble(ID="",
                     GOID=terms[missing],
                     link=ifelse(any(missing),"missing",character(0)))) %>%
    # and collate the ancestors
    bind_rows(
      # and 
      mclapply(
        # 3. function to process the list of obtained offsprings, and collect the 
        # corresponding unique gene ID in the search space
        # this step is parallelised
        sapply(
          # 2. evaluate the expressions
          lapply(
            # 1. construct the expression by indentifying the Ontology (BP,CC or MF)
            # to query the GO database for offsprings
            paste0("as.list(GO",Ontology(terms[!missing]),"OFFSPRING['",terms[!missing],"'])"),
            str2expression),
          eval),
        function(offsprings,associatedIDs){
          associatedIDs[associatedIDs$GOID %in% offsprings,"ID"] %>% unique
        },associatedIDs,mc.cores=mc.cores) %>% 
        # convert the nested list into a nested tibble
        enframe(name=c("GOID")) %>% 
        # unnest (same as melting)
        unnest(cols=c(value)) %>%
        # add an annotation column
        mutate(link="ancestor")
    )
}

