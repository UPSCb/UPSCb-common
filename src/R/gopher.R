require(tidyr)
require(tibble)
require(dplyr)
require(jsonlite)

gopher <- function(genes=character(0),
                   background = NULL,
                   task = list("go",
                            "pfam",
                            "kegg"),
                   alpha=0.05,
                   host="https://microasp.upsc.se",
                   port=5432,
                   url="json",
                   endpoint=c("enrichment","gene-to-term","get-sets","term-to-gene")) {

  # arguments
  if(!is.list(task)){
    task <- list(task)
  }
  
  endpoint <- match.arg(endpoint)
  
  # sanity check
  # any transcript IDs?
  # CAVEAT: very UPSC specific (almost)
  # replace by a function that gets the complete pop from gopher and use that to test
  if(any(grepl("\\.\\d+$",genes))){
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
  
  # request
  request <- httr::POST(paste0(host,":",port,"/",url,"/",endpoint),
             body = list(
               genes = genes,
               background = background,
               target = task,
               alpha = alpha
             ),
             encode = "json")
  
  # implement me
  # request <- httr::POST(paste0("https://microasp.upsc.se",":","5432","/","potra","/","term-to-gene"),
  #                       body = list(
  #                         target=c(
  #                           list(
  #                             name = "go",
  #                             terms= I(c("GO:0015979"))
  #                           )
  #                         )
  #                       ),
  #                       encode = "json")
  # 
  
  # process
  parsed <- jsonlite::fromJSON(httr::content(request, as = "text",
                                             encoding = "UTF-8"))
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
    if(is.character(f)) {
      message(f)
      return(NULL)
    }
    return(as.tibble(f) %>% arrange(padj))
  }, tasks)
  names(parsed) <- tasks
  return(parsed)
}
