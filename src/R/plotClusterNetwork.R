#plot the nodes of a network arguments are Module, counttable, adjacency matrix, label vector
#Requires arguments as follows:
# <mod> is a module number which is checkagainst the label vector
# <frame> the original data frame or matrix with the gene data. Colnames must be present
# <matx> an adjacency matrix. If none is provided it can be calculated using the code in line 17 and might require a LOT of time
# <labels> vector carrying the module labels
# <method> see ?graph.adjaceny

require(igraph)
require(WGCNA)

extractGenes <- function(mod,frame,labels){
  colnames(frame)[labels == mod]
}

checkNode <- function(mod,frame,matx,labels,method='min',threshold=0.5,degrees=1){
  #matx = adjacency(matx,power=31,type="signed")
  mat = as.matrix(matx[extractGenes(mod,frame,labels),extractGenes(mod,frame,labels)])
  mat <- ifelse(abs(mat)<threshold,0,mat)
  for(i in 1:nrow(mat)){
    mat[i,i]<-0
  }

  graph = graph.adjacency(mat,mode=method,weighted=TRUE,)
  graph = delete.vertices(graph,degree(graph)<degrees)

  print(paste(length(V(graph)),'of','(', round(length(V(graph))/length(extractGenes(mod,frame,labels)),4)*100,'% )' ,length(extractGenes(mod,frame,labels)),'vertices with',degrees,'or more degrees and at least', threshold, 'adjacency'))

  v = V(graph)$name

  lab = TRUE
  maxWeight = 2
  siz = 7

  if(length(v)>50){
    lab = FALSE
    siz = 5
    maxWeight = 2
  }

  if(lab==FALSE){igraph.options(vertex.label=c(1:length(v)),vertex.label.cex=0.75)} else {igraph.options(vertex.label=v)}

  igraph.options(edge.curved=FALSE,vertex.size=siz,edge.width=c(1.5**E(graph)$weight*maxWeight),edge.color=NULL,vertex.label.color='black')
  plot.igraph(graph)
  return(list('vertices' = v))
}
