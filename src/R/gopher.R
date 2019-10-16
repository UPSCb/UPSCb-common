
# Ex. y <- gopher(c("GO:0034641", "GO:0071705"), url="pabies", task=c("go"), endpoint = "term-to-gene")
# y$`GO:0071705`
# [1] "MA_10427660g0010" "MA_10436496g0010" "MA_552029g0010"   "MA_59139g0010"    "MA_847283g0010"  
# y$`GO:0034641`
#[1] "MA_10053307g0010" "MA_10263024g0010" "MA_10426126g0010" "MA_15844g0010"    "MA_169451g0010"   "MA_3184966g0010"  "MA_40229g0010"   

require(tidyr)
require(tibble)
require(dplyr)
require(jsonlite)

gopher <- function(genes=character(0),
                   background = NULL,
                   task = list("go",
                            "pfam",
                            "kegg",
                            "mapman"),
                   alpha=0.05,
                   host="https://watson.plantphys.umu.se",
                   port=5432,
                   url="json",
                   endpoint=c("enrichment","gene-to-term","get-sets","term-to-gene")) {

  
  # arguments
  if(endpoint != "term-to-gene" && !is.list(task)){
    task <- as.list(task)
  }
  
  endpoint <- match.arg(endpoint)
  
  #term-to-gene behavior
  if(endpoint=="term-to-gene") {
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

  if(endpoint=="term-to-gene"){
    body = list(target = list(terms = enrTerms, name = task))
  } else {
    body = list(
      genes = genes,
      background = background,
      target = task,
      alpha = alpha
    )
  }
  
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
  if(is.null(request))
    return(NULL)

  
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
  
  tasks <- names(parsed)
  # return
  parsed <- lapply(seq_along(parsed), function(i, task){
    f <- parsed[[i]]
    n <- task[i]
    if(length(f) == 0 | is.null(f))
    {
      message("No enrichments found in task: ", n)
      return(NULL)
    }
    if(! is.null(f$err)) {
      message(f$err)
      return(NULL)
    }
    return(as_tibble(f) %>% arrange(padj))
  }, tasks)
  names(parsed) <- tasks
  return(parsed)
}
