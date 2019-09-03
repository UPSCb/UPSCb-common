library(genomeIntervals)
library(GenomicRanges)

"detailNonModelRanges" <- function (
  incoming=incoming, 
  gff3=gff3, 
  #orgdb,
  dist = 5000,
  promoter = c(3000, 1000),
  max.intron = 1.1e+05
  #key.field = "ENTREZID"
  #name.field = "SYMBOL" 
  ){
  
  ## read the gff3
  gff3 <- readGff3(gff3)
  
  ## convert to GRanges
  curex <- as(gff3[gff3$type=="exon",],"GRanges")
  
  ## extract some info
  gene.id <- as.character(curex$Parent)
  gene.str <- as.logical(strand(curex) == "+")
  
  ## and remove the metadata
  elementMetadata(curex) <- NULL
  
  ## Ignore the annotation
  #anno <- select(orgdb, keys = gene.id, columns = name.field, 
  #               keytype = key.field)
  
  #anno
  
  ## initialise object and record the genes
  ## this could also be done using a factor
  n.entries <- length(gene.id)
  
  #if (nrow(anno) != n.entries) {
  #  stop("possible many-to-one relationship between key and name fields")
  #}
  #all.names <- list()
  #do.check <- !key.field %in% name.field
  #for (x in 1:length(name.field)) {
  #  cur.name <- anno[[name.field[x]]]
  #  if (!is.character(cur.name)) {
  #    cur.name <- as.character(cur.name)
  #  }
  #  if (do.check) {
  #    cur.name <- ifelse(is.na(cur.name), paste0("<", anno[[key.field]], 
  #                                               ">"), cur.name)
  #  }
  #  all.names[[x]] <- cur.name
  #}
  #gene.name <- do.call(paste, c(all.names, sep = ";"))
  
  gene.name <- gene.id
  
  is.diff <- c(TRUE, gene.id[-1] != gene.id[-n.entries] | diff(as.integer(seqnames(curex))) != 
                 0L | diff(gene.str) != 0L | start(curex)[-1] - end(curex)[-n.entries] - 
                 1L > max.intron)
  
  ## collate the exon data
  gene.id <- cumsum(is.diff)
  ngenes <- sum(is.diff)
  output <- .Call(csaw:::cxx_collate_exon_data, gene.id, gene.str, 
                  start(curex), end(curex))
  if (is.character(output)) {
    stop(output)
  }
  ex.num <- output[[1]]
  gb.collected <- output[[2]]
  gb.ref <- gb.collected[[1]]
  gb.ranges <- GRanges(seqnames(curex)[gb.ref], 
                       IRanges(gb.collected[[2]], 
                               gb.collected[[3]]), 
                       strand = !gene.str[gb.ref], seqinfo = seqinfo(curex))
  names(gb.ranges) <- names(curex)[gb.ref]
  if (length(promoter) != 2L) {
    stop("need upstream/downstream specification in promoter argument")
  }
  prom.ranges <- suppressWarnings(trim(promoters(gb.ranges, 
                                                 upstream = promoter[1], downstream = promoter[2])))
  names(prom.ranges) <- names(gb.ranges)
  curex <- c(curex, prom.ranges, gb.ranges)
  gene.name <- c(gene.name, gene.name[gb.ref], gene.name[gb.ref])
  gene.id <- c(gene.id, gene.id[gb.ref], gene.id[gb.ref])
  ex.num <- c(ex.num, integer(length(gb.ref)), rep(-1L, length(gb.ref)))
  gene.str <- c(gene.str, gene.str[gb.ref], gene.str[gb.ref])
  if (missing(incoming)) {
    curex$symbol <- gene.name
    curex$exon <- ex.num
    curex$internal <- gene.id
    curex
    stop("exit here")
  }
  if (any(strand(incoming) != "*")) {
    warning("strandedness in incoming regions is ignored when overlapping")
    strand(incoming) <- "*"
  }
  full.lap <- findOverlaps(incoming, curex)
  flank.only <- ex.num > 0L
  to.flank <- curex[flank.only]
  left.flank <- suppressWarnings(trim(flank(incoming, dist)))
  left.lap <- findOverlaps(left.flank, to.flank)
  left.dist <- start(incoming)[queryHits(left.lap)] - end(to.flank)[subjectHits(left.lap)]
  left.nolap <- left.dist > 0L
  right.flank <- suppressWarnings(trim(flank(incoming, dist, 
                                             start = FALSE)))
  right.lap <- findOverlaps(right.flank, to.flank)
  right.dist <- start(to.flank)[subjectHits(right.lap)] - end(incoming)[queryHits(right.lap)]
  right.nolap <- right.dist > 0L
  all.strs <- .Call(csaw:::cxx_annotate_overlaps, length(incoming), 
                    queryHits(full.lap) - 1L, subjectHits(full.lap) - 1L, 
                    queryHits(left.lap)[left.nolap] - 1L, which(flank.only)[subjectHits(left.lap)][left.nolap] - 
                      1L, left.dist[left.nolap], queryHits(right.lap)[right.nolap] - 
                      1L, which(flank.only)[subjectHits(right.lap)][right.nolap] - 
                      1L, right.dist[right.nolap], gene.name, ex.num, gene.id, 
                    gene.str)
  if (is.character(all.strs)) {
    stop(all.strs)
  }
  names(all.strs) <- c("overlap", "left", "right")
  return(all.strs)
}
