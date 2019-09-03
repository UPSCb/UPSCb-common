## libs
library(genomeIntervals)
library(IRanges)

## wd
setwd("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/gff")
## setwd("/home/paulina/20_ColdAcclimation_Project/")

## read the GFF
gff <- readGff3("Ptrichocarpa_v3.0_210_gene_exons.gff3")
## gff <- readGff3("/home/paulina/references/arabidopsis/TAIR10_GFF3_genes.gff")

## preprocess the file
## exonAttrs <- annotation(gff)$gffAttributes[gff$type == "exon"]
# annotation(gff)$gffAttributes[gff$type == "exon"] <- paste(paste(sub("Parent","ID",exonAttrs),
#                                                                  "exon",
#                                                                  unlist(lapply(runLength(Rle(sub("Parent=","",exonAttrs))),":",1),
#                                                                         use.names=FALSE),sep="."),exonAttrs,sep=";")

na.omit(gff[getGffAttribute(gff,"Parent")=="AT5G58190.1",])

## get the gene <-> pacid map
sel <- gff$type == "mRNA"
idMap <- data.frame(getGffAttribute(gff[sel],"ID"),
                    getGffAttribute(gff[sel],"Parent"))

## extract the exons and group by gene ID
sel <- gff$type == "exon"
mRnaID <- getGffAttribute(gff[sel],"Parent")

## avoid tRNA, rRNA, etc. 
rngs <- IRanges(start=gff[sel,1],end=gff[sel,2])[mRnaID %in% idMap$ID]

## create a set of synthetic exons
rngList <- IRanges::reduce(split(rngs,
                                 idMap[match(mRnaID[mRnaID %in% idMap$ID],idMap$ID),"Parent"]))

## export the gene, mrna and exon as a gff3
## create the new gff object
## select the gene
sel <- gff$type == "gene"

## create the gene gff
geneID <- getGffAttribute(gff[sel],"ID")
geneGff <- gff[sel][geneID %in% idMap$Parent]

## create the mRNA gff
mRnaGff <- gff[sel][geneID %in% idMap$Parent]
mRnaGff$type <- "mRNA"
mRnaGff$gffAttributes <- paste(paste(sub(";",".0;",mRnaGff$gffAttributes),"0;Parent=",sep="."),
                               geneID[geneID %in% idMap$Parent],sep="")

## create the exon gff
rngList <- rngList[match(geneID[geneID %in% idMap$Parent],names(rngList))]
exonNumber <- elementLengths(rngList)
exonGff <- gff[rep(which(sel)[geneID %in% idMap$Parent],exonNumber)]
exonGff[,1] <- unlist(start(rngList))
exonGff[,2] <- unlist(end(rngList))

exonID <- sapply(exonNumber,":",1)
sel <- geneGff$strand == "+"
exonID[sel] <- sapply(exonID[sel],rev)
ID <- getGffAttribute(exonGff,"ID")
exonGff$gffAttributes <- paste("ID=",paste(ID,"exon",unlist(exonID,use.names=FALSE),sep="."),
                               ";Name=",paste(ID,"exon",unlist(exonID,use.names=FALSE),sep="."),
                               ";Parent=",paste(ID,"0",sep="."),sep="")
exonGff$type <- "exon"

## combine
newgff <- c(geneGff,mRnaGff,exonGff)

## change the source
newgff$source <- "UPSC"

## sort
newgff <- newgff[order(seq_name(newgff),newgff[,1],factor(as.character(newgff$type),labels=c(1:3),levels=c("gene","mRNA","exon"))),]

## write
writeGff3(newgff,file="Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3")
## writeGff3(newgff,file="/home/paulina/references/arabidopsis/TAIR10_synthetic-transcripts.gff3")


## done
