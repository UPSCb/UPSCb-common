#' ---
#' title: "Zygotic Embryogenesis first degree neighbour"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(igraph)
})

#' * Data
edgelist <- read.table(here("share/Day5/seidr/backbone-edgelist.txt"),
                       header=FALSE,as.is=TRUE)

#' # Network
#' We create a graf
graf <- graph.edgelist(as.matrix(edgelist))

#' We can check the structure of the graf
#' 
#' There are 3 subgrafs, one with almost all genes, and 2 involving only 2 genes each
#' (look at the $no and $csize list elements)
clusters(graf)

#' From the differential expression results, we can load our DE genes
#' 
#' Get them from 
#' 
#' http://franklin.upsc.se:3000/materials/Exercises/Day5/geneNetworkPreparation.zip
goi <- sub("\\.1","",read.csv("DE-genes.csv",as.is = TRUE)[,1])

#' Extract the first degree neighbours of these genes from the network
subgraf <- make_ego_graph(graf,1,
                          get.vertex.attribute(graf,"name") %in% goi)

barplot(table(sapply(lapply(subgraf,clusters),"[[","csize")),
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")

#' combine all these networks together
fdn <- Reduce("%u%",subgraf)

#' Look at how many clusters we get and how many nodes are involved
clusters(fdn)

#' Let's export the data for visualisation
write_graph(fdn,format = "graphml",file="firstDegreeNeighbour.graphml")

#' # Session Info
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```

