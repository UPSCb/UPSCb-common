stopifnot(
  suppressPackageStartupMessages({
    require(tidyverse)
    require(Biostrings)
  })
  )

## GC percent per gene from a fasta input with gene wise sequences
## Gene name to be truncated until first space

gc_from_fasta <- function(fasta_file) {
  # check input
  stopifnot(is.character(fasta_file), length(fasta_file) == 1)
  
  # read fasta
  fasta <- Biostrings::readDNAStringSet(fasta_file)
  
  # base frequencies
  freq <- Biostrings::alphabetFrequency(
    fasta,
    baseOnly = TRUE,
    collapse = FALSE
  )
  
  GC_percent <- freq %>% 
    as.data.frame(row.names = sub(" .*", "", names(fasta))) %>% 
    rownames_to_column(var = "Gene") %>% 
    as_tibble() %>% 
    rowwise() %>%
    mutate(GC = sum(G, C),
           AGCT = sum(A, G, C, T),
           GC_percent = (GC/AGCT)*100) %>%
    dplyr::select(Gene, GC_percent) %>% 
    as_data_frame()
  
  return(GC_percent)
}
