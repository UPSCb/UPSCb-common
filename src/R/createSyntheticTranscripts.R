"createSyntheticTranscripts" <- function(gff3,
                                         output = c("Genome_intervals","GRanges"),
                                         features = c("mRNA", "tRNA", "miRNA")) {
  require(genomeIntervals)
  require(S4Vectors)

  features <- match.arg(features, several.ok = TRUE)
  output <- match.arg(output)

  gff <- readZeroLengthFeaturesGff3(gff3)

  ## get the gene <-> pacid map
  # This is mRNA IDs and their parents (genes)
  sel <- gff$type %in% features
  idMap <- data.frame(type = gff[sel]$type,
                      getGffAttribute(gff[sel],"ID"),
                      getGffAttribute(gff[sel],"Parent"))

  ## extract the exons and group by gene ID
  sel <- gff$type == "exon"

  ## we can drop multiple Parents (i.e. comma separated Parent values as we are
  ## collapsing them anyway)
  mRnaID <- sub(",.*","",getGffAttribute(gff[sel],"Parent"))

  ## avoid unwanted features
  rngs <- IRanges::IRanges(start = gff[sel, 1],
                           end = gff[sel, 2])[mRnaID %in% idMap$ID]

  ## create a set of synthetic exons
  rngList <- IRanges::reduce(
    split(rngs, idMap[match(mRnaID[mRnaID %in% idMap$ID], idMap$ID), "Parent"]))

  ## export the gene, exon and features as gff3
  ## create the new gff object
  ## select the gene
  sel <- gff$type == "gene"

  ## create the gene gff
  geneID <- getGffAttribute(gff[sel],"ID")
  geneGff <- gff[sel][geneID %in% idMap$Parent]

  ## create gffs for each feature
  featureGff <- Reduce(c, lapply(features, function(f) {
    fGff <- gff[sel][geneID %in% idMap$Parent[idMap$type == f]]
    fGff$type <- f
    fGff$gffAttributes <- paste(paste(
      sub(";", ".0;", fGff$gffAttributes), "0;Parent=", sep="."),
      geneID[geneID %in% idMap$Parent[idMap$type == f]], sep="")
    fGff
  }))

  ## create the exon gff
  rngList <- rngList[match(geneID[geneID %in% idMap$Parent], names(rngList))]
  exonNumber <- elementNROWS(rngList)
  exonGff <- gff[rep(which(sel)[geneID %in% idMap$Parent], exonNumber)]
  exonGff[,1] <- unlist(start(rngList))
  exonGff[,2] <- unlist(end(rngList))

  exonID <- sapply(exonNumber, ":", 1)
  sel <- geneGff$strand == "+"
  exonID[sel] <- sapply(exonID[sel], rev)
  ID <- getGffAttribute(exonGff, "ID")
  exonGff$gffAttributes <- paste0("ID=", paste(ID, "exon", unlist(exonID, use.names=FALSE), sep="."),
                                 ";Name=", paste(ID, "exon", unlist(exonID, use.names=FALSE), sep="."),
                                 ";Parent=", paste(ID,"0",sep = "."))
  exonGff$type <- "exon"

  ## combine
  newgff <- c(geneGff, featureGff, exonGff)

  ## change the source
  newgff$source <- "UPSC"

  ## sort
  newgff <- newgff[order(seqnames(newgff), newgff[, 1],
                         factor(as.character(newgff$type),
                                labels = seq_len(2 + length(features)),
                                levels = c("gene", features, "exon"))), ]

  return(switch(output,
                "Genome_intervals" = newgff,
                "GRanges" = as(newgff, "GRanges"),
                stop(paste("Cannot generate output of type", output))))
}
