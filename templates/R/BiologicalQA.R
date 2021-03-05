#' ---
#' title: "CHANGEME Biological QA"
#' author: "CHANGEME"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Sample information
#' ```{r CHANGEME1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
samples <- read_csv(here("doc/CHANGE-ME.csv"),
                      col_types=cols(.default=col_factor()))

#' tx2gene translation table
#' ```{r CHANGEME2,eval=FALSE,echo=FALSE}
#' # This file is necessary if your species has more than one transcript per gene.
#' #
#' # It should then contain two columns, tab delimited, the first one with the transcript
#' # IDs and the second one the corresponding
#' #
#' # If your species has only one transcript per gene, e.g. Picea abies v1, then
#' # comment the next line
#' ```
tx2gene <- suppressMessages(read_delim(here("reference/annotation/CHANGE-ME"),delim="\t",
                                 col_names=c("TXID","GENE")))

#' # Raw data
filelist <- list.files(here("data/salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info
#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are inthe same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
stopifnot(all(match(sub("_sortmerna.*","",basename(dirname(filelist))),
                    samples$SampleID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$SampleID

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, replace with the following
#' counts <- suppressMessages(round(tximport(files = filelist, type = "salmon",txOut=TRUE)$counts))
#' ```
txi <- suppressMessages(round(tximport(files = filelist, 
                                  type = "salmon",
                                  tx2gene=tx2gene)))
counts <- txi$counts

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by CHANGEME
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=CHANGEME)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by CHANGEME. 
#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely indentify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(CHANGEME=samples$CHANGEME[match(ind,samples$CHANGEME)]) %>% 
  mutate(SamplingTime=samples$SamplingTime[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=CHANGEME)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ CHANGEME)

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Assess whether there might be a difference in library size linked to a
#' given metadata
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This would need to be ggplot2'ed
#' ```
boxplot(split(sizes,dds$CHANGEME),las=2,
        main="Sequencing libraries size factor by Tissue")

plot(sizes,log10(colSums(counts(dds))),ylab="log10 raw depth",xlab="scaling factor",
     col=rainbow(n=nlevels(dds$CHANGEME))[as.integer(dds$CHANGEME)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$CHANGEME)),
       legend=levels(dds$CHANGEME),cex=0.6)


#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=2

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables aas needed.
#' ```
nlevel=nlevels(dds$MDay) * nlevels(dds$MGenotype)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  
#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=CHANGEME,shape=CHANGEME,text=CHANGEME)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Sequencing depth
#' Number of genes expressed per condition at different cutoffs
conds <- factor(paste(dds$CHANGEME,dds$CHANGEME))
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' ### Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 2

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="")

#' ### Hierarchical clustering
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 1000, parallel = TRUE)

#' plot the clustering with bp and au
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' bootstrapping results as a table
print(hm.pvclust, digits=3)


#' ## Conclusion
#' CHANGEME
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
