suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
})

RFAM_VERSION <- "14.4"
HTTPS_SERVER <- file.path("ftp://ftp.ebi.ac.uk/pub/databases/Rfam",RFAM_VERSION,"fasta_files")

# Taxonomy <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
# tar -zxf taxdump.tar.gz names.dmp
# sed -i 's:\t::g' names.dmp
# sed -i 's:"::g' names.dmp
# sed -i 's:|$::g' names.dmp

Tax <- read_delim_chunked("/mnt/picea/storage/reference/Taxonomy/20210226/names.dmp",
                          callback=DataFrameCallback$new(function(chunk,pos){chunk %>% filter(Type=="scientific name") %>% select(c("ID","Name"))}),
                          delim="|",
                  col_names=c("ID","Name","Description","Type"),
                  col_types=cols(ID=col_double(),.default=col_character()))


R5S <- readDNAStringSet(file.path(HTTPS_SERVER,"RF00001.fa.gz"))
  
names(R5S) %>% head %>% str_spli

R5.8S <- readDNAStringSet(file.path(HTTPS_SERVER,"RF00002.fa.gz"))
names(R5.8S)


