### ================================================================
## libraries
### ================================================================
library(RSQLite)
library(igraph)

### ================================================================
## working dir
### ================================================================
#setwd("/mnt/picea/storage/reference/Taxonomy/20160516")
# the working dir is the local dir by default

### ================================================================
## Redo the division table
### ================================================================
con <- dbConnect(dbDriver("SQLite"),dbname="taxonomy.sqlite")
edges <- dbGetQuery(con,"SELECT tax_id, parent_id from node;")

## this lists all superkingdom / kingdom
divs <- dbGetQuery(con,"SELECT tid, nam from taxonomy t LEFT JOIN node n ON t.tid == n.tax_id where rank IN ('superkingdom','kingdom') AND class LIKE 'scientific name%';")

## get the graph (taxonomy is a DAG)
## remove the parents so that it breaks into subgraphs
graph <- graph.edgelist(as.matrix(edges[!edges$parent_id %in% edges[edges$tax_id %in% divs[,1],"parent_id"],]))

## process the divs
## add them to the table
## and update the nodes
apply(divs,1,function(ro,mem){
  dbGetQuery(con, paste("INSERT INTO division VALUES (",ro[1],",'",ro[2],"')",sep=""))
  dbGetQuery(con,paste("UPDATE node set div_id=",ro[1],"WHERE tax_id in (",paste(which(mem == mem[as.integer(ro[1])]),collapse=","),")"))
},clusters(graph)$membership)

## a validation
validation.sql <- paste("sum(case when div_id ==",c(0,divs$tid),"then 1 else 0 end)",collapse=",")
sum(dbGetQuery(con,paste("select", validation.sql, "from node;"))) == dbGetQuery(con,"select count() from node;")

## done
dbDisconnect(con)
