#' ---
#' title: "CHANGEME Biological QA"
#' subtitle: "CHANGEME Project title"
#' author: "CHANGEME"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 3
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#'    includes:
#'      before_body: header.html
#'      after_body: footer.html
#'    css: style.css
#' ---
#' 
#' <hr />
#' &nbsp;
#' 
#' # Setup
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(PCAtools)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(RColorBrewer)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(9,"Blues")

#' # Data
#' * Sample information
#' ```{r Instructions1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
samples <- read_csv(here("doc/CHANGE-ME.csv"),
                    col_types=cols(.default=col_factor()))

#' * tx2gene translation table
#' ```{r Instructions2,eval=FALSE,echo=FALSE}
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

#' * Raw data
filelist <- list.files(here("data/salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' * Sanity check to ensure that the data is sorted according to the sample info
#' ```{r Instructions3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are inthe same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
stopifnot(all(match(sub("_CHANGE-ME.*","",basename(dirname(filelist))),
                    samples$SampleID) == 1:nrow(samples)))

#' * add filelist to samples as a new column
samples %<>% mutate(Filenames = filelist)

#' * export full rank samples
write_tsv(samples,here("doc/samples_full_rank.txt"))

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, or if you are conducting QA 
#' in transcript level replace with the following:
#' txi <- suppressMessages(tximport(files = samples$Filenames, type = "salmon",txOut=TRUE))
#' ```
txi <- suppressMessages(tximport(files = samples$Filenames,
                                 type = "salmon",
                                 tx2gene=tx2gene))
counts <- txi$counts
colnames(counts) <- samples$SampleID
#' 
#' <hr />
#' &nbsp;
#' 
#' # Quality Control
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by CHANGEME
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=CHANGEME)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(CHANGEME), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' `r emoji("point_right")` **We observe almost no difference in the raw sequencing depth**
#' 
#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **The cumulative gene coverage is as expected**
#' 
#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
#' 
#' ## Per-sample expression

dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(CHANGEME=samples$CHANGEME[match(ind,samples$CHANGEME)]) %>% 
  mutate(SamplingTime=samples$SamplingTime[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=CHANGEME)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **All samples have the same sequencing depth characteristics and there is no deviation when we look at one or the other variable**
#' 
#' * Export raw expression data
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
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

colnames(dds) <- samples$SampleID

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' ## size factors 
#' (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **There is almost no differences in the libraries' size factors. They are all within +/1 10% of the average library size.**
#' 
#' Assess whether there might be a difference in library size linked to a
#' given metadata
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This would need to be ggplot2'ed
#' ```
boxplot(split(t(normalizationFactors(dds)),dds$CHANGEME),las=2,
        main="Sequencing libraries size factor by CHANGEME",
        outline=FALSE)

#' `r emoji("point_right")` **The scaling factor distribution is (in)dependent of the CHANGEME variable.**

plot(colMeans(normalizationFactors(dds)),
     log10(colSums(counts(dds))),ylab="log10 raw depth",
     xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$CHANGEME))[as.integer(dds$CHANGEME)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$CHANGEME)),
       legend=levels(dds$CHANGEME),cex=0.6)

#' `r emoji("point_right")` **The scaling factor appear linearly proportional to the sequencing depth.**
#' 
#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:
meanSdPlot(vst[rowSums(vst)>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' `r emoji("point_right")` **We can conclude that the variance stabilization worked adequately**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' ### Scree plot
#' 
#' We define the number of variable of the model:
nvar <- "CHANGEME (for a number)"

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables as needed.
#' ```
nlevel <- nlevels(dds$CHANGEME) * nlevels(dds$CHANGEME)

#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' We plot the percentage explained by different components and try to empirically assess whether
#' the observed number of components would be in agreement with our model's assumptions.
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' * the black dotted, annotate lines represent the optimal number of components 
#' reported by the horn and elbow methods.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))

#' `r emoji("point_right")` **The first component explains 40% of the data variance. Both metrics, Horn and Elbow suggest that one or two components are those that are informative. Indeed the slope of the curve is fairly linear past PC3 and that would indicate that the remaining PCs only capture sample specific noise. While this is only empirical, the scree plot support having only few variables of importance in the dataset.**
#'
#'
#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#' * PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=CHANGEME,shape=CHANGEME,text=CHANGEME)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,
       colby = 'CHANGEME',
       colLegendTitle = 'CHANGEME',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' * PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=CHANGEME,shape=CHANGEME,text=CHANGEME)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,x = 'PC1', y = 'PC3',
       colby = 'CHANGEME',
       colLegendTitle = 'CHANGEME',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ```{r subplot, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```
#'
#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby='CHANGEME',shape='CHANGEME'))

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ### Correlation
#' This is a plot showing the correlation between the PC and the model variables. Note that while this might be relevant 
#' for a linear variable, it is less so for categorical variables. Sorting categorical variables in a linear order according to the PCs above might help.
#' 
#' ```{r beffect, echo=FALSE, eval=FALSE
#' # You could do as follows to reorder a categorial variable
#' p$metadata$Beffect <- ifelse(as.integer(dds$Batch)==1,3,as.integer(dds$Batch)-1)
#' ````
#' 
#' 
#' Plotting only the relevant variables.
#' 
#' ```{r validation, echo=FALSE, eval=FALSE}
#' # Plotting all the variables just as a control, the categorical and the transformed ones.
#' suppressWarnings(eigencorplot(p,metavars=c('Beffect','Day',"Batch","Time")))
#' ```
suppressWarnings(eigencorplot(p,metavars=c('CHANGEME','CHANGEME')))

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ### Samples Distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- dds$SampleID
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=pal)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- factor(paste(dds$CHANGEME,dds$CHANGEME))
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' Plotting the number of genes that are expressed (at any level)
do.call(rbind,split(t(nrow(vst) - colSums(vst==0)),samples$CHANGEME)) %>% as.data.frame() %>% 
  rownames_to_column("CHANGEME") %>% pivot_longer(starts_with("V")) %>% 
  ggplot(aes(x=CHANGEME, y=value, fill=CHANGEME)) + geom_dotplot(binaxis = "y", stackdir = "center")

#' `r emoji("point_right")` **CHANGEME**
#' 
#' ## Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#'
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

#' `r emoji("point_right")` **Here a cutoff of CHANGEME is applied**
vst.cutoff <- 2

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#' `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:
#'
#' ```{r CHANGEME9,eval=FALSE,echo=FALSE}
#' Optionally You can decide which variable to select for annotation. 
#' add the following to your pheatmap call (select the variable you want)
#' This will be a colored ribbon on top of the heatmap for better reading
#' annotation_col = samples %>% select(CHANGEME) %>% as.data.frame()
#' 
#' if you want, you can also modify the colors to better match with previous graphs. 
#' e.g. if you want a colored ribbon for Temperature (16,23,27) (adjust this example based on your need and add to pheatmap call)
#' annotation_colors = list(Temperature = c("23" = "#F8766D", "16" = "#00BF7D", "27" = "#00B0F6"))
#' 
#' Alternatively you can use the heatmap.2 version
#' #heatmap.2 version
#' hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")},
#'                 labRow = NA,
#'                 trace = "none",
#'                 labCol = conds,
#'                 col=hpal)
#' ```
#'  
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)

#' `r emoji("point_right")` **CHANGEME**
#'
#' ## Clustering of samples
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This wouldonly works with the gplots heatmap.2, not the pheatmap
#' plot(as.hclust(hm$colDendrogram),xlab="",sub="")
#' ```
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' `r emoji("point_right")` **CHANGEME**
#' 
#' <details><summary>bootstrapping results as a table</summary>
#' ```{r bootstrapping results as a table}
#' print(hm.pvclust, digits=3)
#' ```
#' </details>
#' 
#' ```{tech rep, echo=FALSE, eval=FALSE}
#' # The block of code is meant to combine tech reps - as it is facultative it is commented out
#' # First create a new variable in your sample object called BioID that identifies uniquely technical replicates, so one value for all tech rep of the same bio rep
#' samples$BioID <- CHANGEME
#' # Merging technical replicates
#' txi$counts <- sapply(split.data.frame(t(txi$counts),samples$BioID),colSums)
#' txi$length <- sapply(split.data.frame(t(txi$length),samples$BioID),colMaxs)
#' # Counts are now in alphabetic order, check and reorder if necessary
#' stopifnot(colnames(txi$counts) == samples$BioID)
#' samples <- samples[match(colnames(txi$counts),samples$BioID),]
#' # Recreate the dds
#' dds <- DESeqDataSetFromTximport(
#'   txi=txi,
#'   colData = samples,
#'   design = ~ Tissue)
#'```
#'
#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` **CHANGE-ME**
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
#' </details>
#'   
#' &nbsp;
#' 
