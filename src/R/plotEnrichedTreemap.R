# Function originally from https://github.com/loalon/Rtoolbox

clusterTreemapColors <- rep(c("#9E1F63","#662D91","#2B3990","#1B75BC","#27AAE1",
                              "#2AB592","#035E31","#009444","#65BC46","#A5CE42",
                              "#F9ED32","#FBB040","#F15A29","#EF4136","#BE1E2D"),15)
clusterTreemapText <- rep(c("white","white","white","white","white",
                            "white","white","white","white","white",
                            "black","black","white","white","white"),15)

plotEnrichedTreemap <- function(x, enrichment = c('go','mapman', 'kegg', 'pfam', 'ko_pathway', 'ko', 'kog','cog'), 
                                namespace = c('none', 'BP', 'MF', 'CC'), 
                                title = "", 
                                de = c("none", "up", "down"), 
                                clusterColor = "#9E1F63", clusterText='black',
                                nameCol="name",
                                namespaceCol="namespace",
                                sizeCol="padj",
                                colorCol="padj", 
                                convertSize=TRUE,legend=TRUE) {
  
  require(treemap)
  require(RColorBrewer)
  
  enrichment <- match.arg(enrichment)
  namespace <- match.arg(namespace)
  de <- match.arg(de)
  
  enrData <- if (is.data.frame(x)){
    x 
  } else if (is.list(x)) {
    x[[enrichment]]
  } else {
    stop("x object must be either a list or a data.frame")
  }
  
  # calculate the size based on padj
  
  enrData$size <- if(convertSize) {
    abs(log10(enrData[[colorCol]]))
  } else {
    enrData[[colorCol]]
  }
  
  #default treemap
  index = nameCol
  fontcolor.labels=clusterText
  fontface.labels=c(2)
  fontsize.labels=c(12)
  inflate.labels = TRUE
  align.labels=list(c("center", "center"))
  position.legend <- if (!legend) {"none"} else {"bottom"}
  #vColor = 'name'
  border.col='black'
  border.lwds=c(7,2)
  palette <- colorRampPalette(c("white", clusterColor))
  palette <- palette(10)
  vColor = c("size")
  if (sizeCol=="padj") {
    vSize = vColor
  } else {
    vSize = sizeCol
  }
  type = "value"
  title.legend="abs(log10(pAdj))"
  bg.labels= 0
  
  if(enrichment =='go' ){
    if(namespace=='none') {
      index = c('namespace', 'name')
      palette = "Set1"
      
      fontcolor.labels=c("#FF000000", "black")
      fontface.labels=c(1, 2)
      fontsize.labels=c(1, 20)
      inflate.labels=TRUE
      align.labels=list(c("left", "top"), c("center", "center") )
      type = "index"
      border.col=c('black', 'white')
      title.legend="GO namespace"
    } else {
      enrData <- enrData[enrData[[namespaceCol]]==namespace,]
    }
  }
  
  # mapman name fix for better visualization
  if(enrichment =='mapman') {
    enrData[[nameCol]] <- gsub("\\."," ",enrData[[nameCol]])
  } 
  
  # Paint it properly if it is up or down regulated
  if(de !='none') {
    if(de=='up') {
      palette = "OrRd"
    }
    if(de=='down') {
      palette = "GnBu"
    }
  } 
  
  # generate treemap
  treemap(enrData, 
          index = index,
          vSize = vSize, 
          palette = palette,
          type = type, 
          vColor = vColor,
          title=title, 
          fontcolor.labels=fontcolor.labels, 
          fontface.labels=fontface.labels,     
          #fontsize.labels=fontsize.labels,
          bg.labels=bg.labels, 
          inflate.labels = inflate.labels ,
          lowerbound.cex.labels = 0, 
          position.legend = position.legend,
          border.col=border.col,         
          border.lwds=border.lwds,
          align.labels=align.labels,
          title.legend=title.legend,
          overlap.labels = 1
  )
}
