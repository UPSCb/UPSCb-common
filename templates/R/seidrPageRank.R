library(here)
library(igraph)
library(readr)

sf <- read_tsv(here("data/seidr/backbone/backbone-2-percent-filtered.txt"),
               col_names=FALSE,col_types=cols(.default=col_character()),
               show_col_types=FALSE)

d.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Directed",1:2]),directed=TRUE)
u1.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Undirected",1:2]),directed=TRUE)
u2.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Undirected",2:1]),directed=TRUE)
graf <- union(d.graf,u1.graf,u2.graf)

pr <- page_rank(graf)$vector
pr
