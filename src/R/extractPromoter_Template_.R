## libs
library(genomeIntervals)
library(Biostrings)

## working directory
setwd("/mnt/picea/storage/reference/Picea-abies/")

## read the gff
gff <- readGff3("GenePrediction_201309/Eugene.gff3")

## read the fa
fa <- readDNAStringSet("fasta/Pabies1.0-genome.fa")

## subset for the sequence present - reduce the mem footprint!
fa <- fa[levels(seq_name(gff))]
gc()

## identify the CDS starting exon and its first codon
## get the CDS annot by strand
p.strand.cds <- gff[gff$type=="CDS" & gff$strand=="+",]
m.strand.cds <- gff[gff$type=="CDS" & gff$strand=="-",]

## get the CDS ID
p.strand.id <- sub("\\.CDS\\.[0-9]+$","",getGffAttribute(p.strand.cds,"ID"))
m.strand.id <- sub("\\.CDS\\.[0-9]+$","",getGffAttribute(m.strand.cds,"ID"))

## find the first CDS of every gene (first on plus strand, last on minus strand)
p.first.cds <- p.strand.cds[match(unique(p.strand.id),p.strand.id),]
m.first.cds <- m.strand.cds[length(m.strand.id) - match(unique(m.strand.id),rev(m.strand.id)) + 1,]

## keep only those that start with an ATG
p.first.cds <- p.first.cds[subseq(fa[seq_name(p.first.cds)],start=p.first.cds[,1],width=3) == DNAStringSet("ATG"),]
m.first.cds <- m.first.cds[reverseComplement(subseq(fa[seq_name(m.first.cds)],end=m.first.cds[,2],width=3)) == DNAStringSet("ATG"),]

## get the promoter ranges and correct for "out of range" cases
p.rngs <- promoters(IRanges(p.first.cds[,1],width=1),downstream=0)
start(p.rngs)[start(p.rngs)<=0] <- 1
m.rngs <- shift(promoters(IRanges(m.first.cds[,2],width=1),upstream=0,downstream=2000),shift=1)
sel <- end(m.rngs) > width(fa[seq_name(m.first.cds)])
end(m.rngs)[sel] <- width(fa[seq_name(m.first.cds)])[sel]

## get the promote sequence
p.prom <- subseq(fa[seq_name(p.first.cds)],start(p.rngs),end(p.rngs))
names(p.prom) <- sub("\\.CDS\\.[0-9]+$","",getGffAttribute(p.first.cds,"ID"))
m.prom <- reverseComplement(subseq(fa[seq_name(m.first.cds)],start(m.rngs),end(m.rngs)))
names(m.prom) <- sub("\\.CDS\\.[0-9]+$","",getGffAttribute(m.first.cds,"ID"))

## write out
prom <- append(p.prom,m.prom)
writeXStringSet(prom,file="fasta/Pabies1.0-promoter-2kb.fa")

## extract the smaller promoter files from it
## 1kb
prom1kb <- subseq(prom,start=ifelse(width(prom)>1000,width(prom)-1000+1,1),end=width(prom))
writeXStringSet(prom1kb,file="fasta/Pabies1.0-promoter-1kb.fa")

## 500bp
prom500b <- subseq(prom,start=ifelse(width(prom)>500,width(prom)-500+1,1),end=width(prom))
writeXStringSet(prom500b,file="fasta/Pabies1.0-promoter-500bp.fa")

