require(here)
require(readr)
require(dplyr)
require(GO.db)


# The following function takes a vector of GO_IDs and returns a vector of GO names

Pull_GO_name_from_ID <- function(GO_IDs){
  go <- keys(GO.db, keytype="GOID")
  
  go_terms <- suppressMessages(AnnotationDbi ::select(GO.db, columns=c("GOID","TERM"), keys=GO_IDs, keytype="GOID")) %>% 
    as_tibble() %>% dplyr::select(TERM) %>% dplyr::pull() %>% paste(collapse = ";")
  
  return(go_terms)
}

## This commented out example shows how to use the function to add a column
## with GO_names to a dataframe that already has a column with GO IDs.


# example_set <- read_csv(here("DE_genes_with_PA_and_TAIR_WithEggnogAnnotation.csv")) %>% 
#  rowwise() %>% 
#  mutate(
#    GO_names = ifelse(grepl("GO:", GOs), Pull_GO_name_from_ID(strsplit(GOs, ",")[[1]]), "-")
#  ) %>% 
#  ungroup()

#write_tsv(example_set, here("DE_genes_with_PA_and_TAIR_WithEggnogAnnotationAndGOnames.tsv"))
