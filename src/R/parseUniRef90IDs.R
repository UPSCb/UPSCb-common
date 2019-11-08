library(tidyverse)

file="/mnt/picea/storage/reference/UniRef90/201908/annotation/uniref90.id"

f <- function(c,i){
  str_match(c,">([^ ]+) .*Tax=(.*) TaxID=(\\d+).*")[,2:4]
}

df <- read_lines_chunked(file,callback=DataFrameCallback$new(f),
                         chunk_size=1e6,progress=TRUE)

saveRDS(df,file="/mnt/picea/storage/reference/UniRef90/201908/annotation/uniref90.id.rds")

write_delim(as.data.frame(df),col.names=FALSE,
            path="/mnt/picea/storage/reference/UniRef90/201908/annotation/uniref90_id-table.txt")
