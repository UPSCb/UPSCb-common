## working dir
setwd("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/annotation")

## read the annot
annotationtable<-read.delim(file="Ptrichocarpa_v3.0_210_annotation_info.txt"
                            header=FALSE,stringsAsFactors=FALSE,row.names=1)

## check that the 3rd column is useless and remove it
all(annotationtable$V3 == annotationtable$V4)
annotationtable$V3 <- NULL

## check the number of gene and transcripts
length(unique(annotationtable$V2))
length(unique(annotationtable$V4))

## set the colnames
colnames(annotationtable) <- c("geneID","trxID","PFAM","PANTHER","KOG","EC","RefSeq","GO","At","Symbol","Description")

## check that some isoforms have different annotations
table(sapply(sapply(split(annotationtable$PFAM,annotationtable$geneID),unique),length))
table(sapply(sapply(split(annotationtable$GO,annotationtable$geneID),unique),length))

## combine these information into a single table
annot<-cbind(geneID=unique(annotationtable$geneID),sapply(colnames(annotationtable)[-c(1:2)],function(co,tab){
  sapply(sapply(split(tab[,co],tab$geneID),unique),paste,collapse="|")
  },annotationtable))

## write it out
write.csv(annot,file="Ptrichocarpa_v3.0_210_gene-annotation.csv")
