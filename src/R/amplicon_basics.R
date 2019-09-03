library("ggplot2")
library("grid")
library("scales")
theme_set(theme_bw())

#alphadiversity:
#logShannon 
map_file="/path_to_map_file/map.txt"
bmsd=import_qiime_sample_data(map_file)
#relevel your mapping file and make sure your categories are factors f.e.:
#bmsd$Date <- relevel(bmsd$Date, "5th.June")
#bmsd$Treatment<-factor(bmsd$Treatment,levels=c("ctrl","fert"))

biom_file_json="/path_to_json_file/json.biom"
biomot=import_biom(biom_file_json,parseFunction=parse_taxonomy_greengenes)
phyloseq=merge_phyloseq(biomot,bmsd)

meth=("Shannon")

erDF <- estimate_richness(phyloseq, split = TRUE, measures = meth)
df <- data.frame(erDF, sample_data(phyloseq))
#check normal distribution of the data - if not f.e. log-transform
df[,paste0("log", meth)] <- log(df[,meth])

#boxplots for alphadiv
theme_set(theme_bw())
ggplot(df, aes(#yourcategory#, logShannon, fill=(#yourcategory#), dodge=#yourcategory#))+stat_boxplot(geom ='errorbar')+
  geom_boxplot()+scale_fill_manual( values = c("blue", "yellow","red"))+ylim(0, 1.5)

#basic statistics
m<-with(df[c(1:12),], mean(logShannon) ) 
sd<-with(df[c(1:12),], sd(logShannon) ) 
se<-sd/sqrt(12)

#top10 genera relative abundance

genera = tax_glom(phyloseq, "Genus")
fs<-transform_sample_counts(genera, function(OTU) OTU/sum(OTU) *100)
topsc <- names(sort(taxa_sums(fs), TRUE)[1:10])
exsc    <- prune_taxa(topsc, fs)
OTU<-otu_table(exsc)
TAX<-tax_table(exsc)
OTU<-as.data.frame(OTU)
OTU_2<-OTU[order(rowSums(OTU),decreasing=T),]
TAX<-as.data.frame(TAX)
TAX_2<-TAX[order(rowSums(OTU),decreasing=T),]
Means<-rowMeans(OTU_2)
SD<-apply(as.matrix(OTU_2),1,sd)
SE<-SD/sqrt(length(OTU_2))
OTU_2$Means<-as.numeric(Means)
OTU_2$SD<-as.factor(SD)
OTU_2$SE<-as.factor(SE)
EMF<-cbind(OTU_2,TAX_2)
write.table(EMF,file="name.csv",sep=",")

#observed alpha box-plot
biom_file_json="/path_to_rarefied_file/json.biom"
biomot=import_biom(biom_file_json,parseFunction=parse_taxonomy_greengenes)
phyloseq_o=merge_phyloseq(biomot,bmsd)
meth=("Observed")

erDFO <- estimate_richness(phyloseq_o, split = TRUE, measures = meth)
dfo <- data.frame(erDFO, sample_data(phyloseq_o))
ggplot(dfo, aes(#yourcategory#, Observed, fill=(#yourcategory#), dodge=#yourcategory#))+stat_boxplot(geom ='errorbar')+
  geom_boxplot()

#Venn-observed
suppressPackageStartupMessages(library(VennDiagram))
Venn<-otu_table(phyloseq_o)
Venn<-as.data.frame(Venn)
s<-Venn[,1:36]
r<-Venn[,37:72]
n<-Venn[,73:108]
s<-s[!!rowSums(s)>=1, ]
r<-r[!!rowSums(r)>=1, ]
n<-n[!!rowSums(n)>=1, ]
rownames(s)->s
rownames(r)->r
rownames(n)->n
plot.new()
#all treatments 
grid.draw(venn.diagram(list(s,
  r,n),scaled=TRUE,
  filename=NULL,lty="blank",
  fill=c("black","brown","green"),category.names=c("s","r","n")))


#PCoA with vst(DESeq2)
data2<-otu_table(phyloseq)

data2<-as.data.frame(data2)
TAX<-tax_table(rootsoil)

conditions <-colnames(data2)
dds <- DESeqDataSetFromMatrix(
  countData = data2,
  colData = data.frame(condition=conditions),
  design = ~ condition)

dds <- estimateSizeFactors(dds)

colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(data2)
vst <- vst - min(vst)

vst<-as.matrix(vst)

OTUvst = otu_table(vst, taxa_are_rows = TRUE)

physeqvst = phyloseq(OTUvst,TAX)
sampledata=sample_data(phyloseq)
physeqvst=merge_phyloseq(physeqvst,sampledata)
ordvst<-ordinate(physeqvst,"PCoA","bray")
p = plot_ordination(physeqvst,ordvst,type="samples", color="yourcategory", shape = "yourcategory")
p=p+geom_point(size=7, alpha=0.75)+scale_color_manual(values = c("blue", "yellow","red"))
p

p2 = plot_ordination(physeqvst, ordvst, type = "split", color = "Phylum", shape = "yourcategory")
p2

