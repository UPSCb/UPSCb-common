## libs
library(genomeIntervals)
library(Biostrings)

## working dir
setwd("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/")

## read the gff
gff <- readGff3(file="gff/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3")
gff <- gff[gff$type=="exon",]

## read the genome
fa <- readDNAStringSet("fasta/Ptrichocarpa_v3.0_210.fa")

## create the gene-model sequences
seqs <- subseq(fa[seq_name(gff)],gff[,1],gff[,2])
names(seqs) <- getGffAttribute(gff,"Parent")

## the Reduce takes for ever (the outter one), try to just use a List
seqs <- Reduce(append,mclapply(unique(names(seqs)),function(nam,seqs){
  DNAStringSet(Reduce("c",seqs[names(seqs) == nam]))
},seqs,mc.cores=16))

## i.e. this might be faster; try out
seqs <- DNAStringSet(mclapply(unique(names(seqs)),function(nam,seqs){
  Reduce("c",seqs[names(seqs) == nam])
},seqs,mc.cores=16))

names(seqs) <- unique(getGffAttribute(gff,"Parent"))
m.sel <- gff[match(unique(names(seqs)),getGffAttribute(gff,"Parent"))]$strand=="-"
seqs[m.sel] <- reverseComplement(seqs[m.sel])

## export them
writeXStringSet(seqs,"fasta/Ptrichocarpa_v3.0_210_synthetic-gene-models.fa")
