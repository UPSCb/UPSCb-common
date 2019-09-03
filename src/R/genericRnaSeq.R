## Nico: That looks nice to me. What you created is an S3 function;
## at one point we might want to change it into an S4 function
## (see http://cran.r-project.org/doc/contrib/Genolini-S4tutorialV0-5en.pdf)
## but that's not the most urgent. Read the PDF whenever you have time,
## it is an interesting read and is good for learning more about
## programming in general


#rm(list=ls())
#rm(samples)

#quit()

## 4) We should now try to convert your script into an executable script
## that means a script that can read arguments. Have a look there to get
## inspiration: rePairFastq.R in the same dir as this file



# samples_csv <- "samplesTest.csv"
# output_prefix <- "Prefix"
# wd_template <- "~/Git/UPSCb/"
# HTSeq_dir <- "~/Git/UPSCb/HTSeq"
# HTSeq_pattern <- "*.txt"
# HTSeq_gsub <- "_.*_STAR\\.txt"
# HTSeq_mc_cores <- 9
# reads_1_or_2 <-2 #or 1 
# pca1<-1
# pca2<-2
# pca3<-3
# Metrics <-"Place"

#rm(list=ls())
#install.packages("optparse", repos="http://R-Forge.R-project.org")
suppressPackageStartupMessages(library(optparse))

#install.packages("shiny")
suppressPackageStartupMessages(library(shiny))


### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(vsn))#Old model? vsn2
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vcd))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(arrayQualityMetrics))

source("~/Git/UPSCb/src/R/plot.multidensity.R")




genericRnaSeq <- function(output_prefix="", samples_csv = "samples.csv",
                       wd_template = "~/Git/UPSCb",HTSeq_dir = "~/Git/UPSCb/HTSeq",
                       HTSeq_pattern = "*.txt",HTSeq_gsub = "_.*_STAR\\.txt",
                       HTSeq_mc_cores = 9,reads_1_or_2 = 2,pca1 = 1,pca2 = 2,
                       pca3 = 3,Metrics = "Place", pal=brewer.pal(8,"Dark2"),
                       norm.fun = "vst")
{

  
 ### ==============================
 ## set the working dir
 ### ==============================
 setwd(dirname(HTSeq_dir))
 getwd()
 
 ### ==============================
 ## Checking the samples file for existence and number of columns
 ### ==============================
 message("Checking the samples file")
 smp <- file.exists(samples_csv)
 stopifnot(smp)
 if (Biobase::testBioCConnection() == TRUE){
   setwd(wd_template)
   comp1 = scan("Template.csv", what = character(), sep = ",", nmax = 6)
 }else{
   comp1 <- c("ExperimentTitle,SampleName,SampleDescription,SequencingDate,FileName,FileLocation")
 }
 comp2 <- unlist(comp1)
 scan1 <- scan(samples_csv, what = character(), sep = ",", nmax = 6)
 scan2 <- unlist(scan1) 
 stopifnot(all(scan2 %in% comp2))
 
 ### ==============================
 ## read the HTSeq files in a matrix 
 ### ==============================
 message("Reading the HTSeq files in a matrix")
 res <- mclapply(dir(HTSeq_dir,pattern=HTSeq_pattern,full.names=TRUE),function(fil){
   read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
 },mc.cores=HTSeq_mc_cores)
 names(res) <- gsub(HTSeq_gsub,"",dir(HTSeq_dir,pattern=HTSeq_pattern))
 res
 

 
 ### ==============================
 ## get the count table
 ### ==============================
 message("Getting the count table")
 addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
 sel <- match(addInfo,res[[1]][,1])
 count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
 colnames(count.table) <- names(res)
 rownames(count.table) <- res[[1]][,1][-sel]
 write.csv(count.table,file=paste(output_prefix,"countTable.csv",sep=""))
 
 
 ### ==============================
 ## extract the HTSeq stat lines  
 ### ==============================
 message("Extracting the HTSeq stat lines")
 count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
 colnames(count.stats) <- names(res)
 rownames(count.stats) <- addInfo # sub("__","",addInfo) sometimes substitution has been performed at this step
 count.stats <- rbind(count.stats,aligned=colSums(count.table))
 count.stats <- count.stats[rowSums(count.stats) > 0,]
 
 ## as percentages
 apply(count.stats,2,function(co){round(co*100/sum(co))})
 
 ### ==============================
 ## plot the stats
 ### ==============================
 pal=brewer.pal(10,"Paired")[1:nrow(count.stats)]
 par(mar=c(3.1,2.1,3.1,0.1)) #par(mar=c(7.1,5.1,4.1,2.1))
 barplot(as.matrix(count.stats),col=pal,beside=TRUE,las=2,main="read proportion",
         ylim=range(count.stats)+ c(0,2e+6))
 legend("topleft","left",fill=pal,legend=rownames(count.stats),bty="n",cex=0.5)
 
 # ### ==============================
 # ## % of the genes are not expressed
 # ## out of a total of NR
 # ### ==============================
 sel <- rowSums(count.table) == 0
 message(sprintf("%s percent",round(sum(sel) * 100/ nrow(count.table),digits=1)))
 message(sprintf("of %s genes are not expressed",nrow(count.table)))
 

 ### ==============================
 ## display the per-gene mean expression
 ## i.e. the mean raw count of every 
 ## gene across samples is calculated
 ## and displayed on a log10 scale
 ### ==============================
 message("Display the per-gene mean expression")
 plot(density(log10(rowMeans(count.table))),col=pal[1],
      main="mean raw counts distribution",
      xlab="mean raw counts (log10)")
 
 

 
 
 
 ### ==============================
 ## The same is done for the individual
 ## samples colored by sample type
 ### ==============================
 pal=pal
 plot.multidensity(log10(count.table),col=sample(pal, size=length(res), replace=TRUE),
                   legend.x="topright",legend.cex=0.5,
                   main="sample raw counts distribution",
                   xlab="per gene raw counts (log10)",lwd=3)
 
 
 
 
 
 
 ### ==============================
 ## For visualization, the data is
 ## submitted to a variance stabilization
 ## transformation using DESeq2. The 
 ## dispersion is estimated independently
 ## of the sample type
 ### ==============================
 samples <- read.csv(samples_csv)
 conditions <- names(res) 
 dds <- DESeqDataSetFromMatrix(countData = count.table,
                               colData = data.frame(condition=conditions),
                               design = ~ condition)
 colData(dds)$condition <- factor(colData(dds)$condition,
                                  levels=unique(conditions))
 vsd <- varianceStabilizingTransformation(dds, blind=TRUE) 
 
 
 ## create the dds object
 dds <- DESeqDataSetFromMatrix(
   countData = count.table,
   colData = data.frame(condition=conditions
   ),
   design = ~ condition)
 
 
 

 ### ==============================
 ## plot the mean against the variance (sd)
 ## it should be flat if the VST worked properly
 ## it is not but it is likely due to the 
 ## variability of the gene expression
 ## This is not a concern for visualization
 ## but might be for normalization
 ### ==============================
 message("Plot the mean against the variance")
 meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))
 title(main="Mean vs Variance - no cutoff")
 
 
 ### ==============================
 ## check the size factors
 ## no big variation, can go for vsd
 ## over rld
 ### ==============================

 message("Check the size factors")
 dds <- estimateSizeFactors(dds)
 sizes <- sizeFactors(dds)
 names(sizes) <- names(res)
 sizes
 
 
 if (norm.fun == "vst"){## do the vst
   message("The vst is computed")
   vst <- assay(vsd)
   colnames(vst) <- colnames(count.table)
   vst <- vst - min(vst)
   write.csv(vst,paste(output_prefix,"vst.csv",sep=""))
   
   
   ## perform a Principal Component
   ## Analysis (PCA) of the data
   ## to do a quick quality assessment
   ## i.e. replicate should cluster
   ## and the first dimensions should 
   ## be explainabled by biological means.
   ### ==============================
   message("PCA")
   pc <- prcomp(t(assay(vsd)))
   percent <- round(summary(pc)$importance[2,]*100)
   smpls <- unique(conditions)
   cf <- substr(x=colnames(count.table),5,5)
   
  
   ### ==============================
   ## Finding the variable list for the PCA plots
   ### ==============================
   
   coln <- colnames(samples)[7:ncol(samples)]  
   col_finder <- function(col1){
     length(unique(samples[,col1]))*2 < length(samples[,col1])/reads_1_or_2
   }
   colsAll <- sapply(coln,col_finder)
   var.list <- coln[grep(T,colsAll)]
   
  
   ### ==============================
   ## Plot the first two  dims
   ### ==============================
   plot(pc$x[,pca1],
        pc$x[,pca2],
        xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
        ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
        col=1,pch=1,
        main="Principal Component Analysis",sub="variance stabilized normalized counts")
   
   
   
   
   ### ==============================
   ##  Plot 3 choosed dimensions of the PCA 
   ### ==============================
   
   
   pca.plot <- function(var,var2){
     if (length(unique(samples[,var])) < 26){
       if (reads_1_or_2 == 2){
         
         pch = c(1:length(unique(samples[,var])))[as.integer(as.factor(samples[,var]))][c(T,F)]
         colours = c(1:length(unique(samples[,var2])))[as.integer(as.factor(samples[,var2]))][c(T,F)]
       }else{
         pch = c(1:length(unique(samples[,var])))[as.integer(as.factor(samples[,var]))]
         colours = c(1:length(unique(samples[,var2])))[as.integer(as.factor(samples[,var2]))]
       }
     }
     
     
     scatterplot3d(pc$x[,pca1],
                   pc$x[,pca2],
                   pc$x[,pca3],
                   xlab=paste("Comp. ",pca1," (",percent[pca1],"%)",sep=""),
                   ylab=paste("Comp. ",pca2," (",percent[pca2],"%)",sep=""),
                   zlab=paste("Comp. ",pca3," (",percent[pca3],"%)",sep=""),
                   color=colours, pch= pch) #pchs <- c(17,19)[as.integer(sex)])
     legend(xy.coords(-3.5,7),pch= pch,col=colours,legend=smpls,bty="n",cex = 0.5) # rep(c(20,17,length(pc$x[,3])),8)
     mar=c(3.1,2.1,3.1,0.1)
     par(mar=mar)
   }
   
   out<-combn(var.list,2)
   mapply(pca.plot,out[1,],out[2,],USE.NAMES=FALSE)
 }else if(norm.fun=="rld"){
   
   ## do the rld
   message("The rld is computed")
   rld <- rlogTransformation(dds)
   rlt <- assay(rld)
   colnames(rlt) <- colnames(count.table)
   write.csv(rlt,paste(output_prefix,"rlt.csv",sep=""))
   
   
   
   ### ==============================
   ## Perform a Principal Component Analysis
   ## on the rld data
   ### ==============================
   message("PCA")
   pc.rld <- prcomp(t(assay(rld)))
   percent.rld <- round(summary(pc.rld)$importance[2,]*100)
   
   ### ==============================
   ## plot the first 3 dimensions
   ## The PCA looks good. first dimension
   ## is clearly developmental gradient
   ## that explains 37% of the data
   ### ==============================
   
   pcrld.plot <- function(var,var2){
     if (length(unique(samples[,var])) < 26){
       if (reads_1_or_2 == 2){
         
         pch = c(1:length(unique(samples[,var])))[as.integer(as.factor(samples[,var]))][c(T,F)]
         colours = c(1:length(unique(samples[,var2])))[as.integer(as.factor(samples[,var2]))][c(T,F)]
       }else{
         pch = c(1:length(unique(samples[,var])))[as.integer(as.factor(samples[,var]))]
         colours = c(1:length(unique(samples[,var2])))[as.integer(as.factor(samples[,var2]))]
       }
     }
     
     
     scatterplot3d(pc.rld$x[,pca1],
                   pc.rld$x[,pca2],
                   pc.rld$x[,pca3],
                   xlab=paste("Comp. ",pca1," (",percent[pca1],"%)",sep=""),
                   ylab=paste("Comp. ",pca2," (",percent[pca2],"%)",sep=""),
                   zlab=paste("Comp. ",pca3," (",percent[pca3],"%)",sep=""),
                   color=colours, pch= pch) #pchs <- c(17,19)[as.integer(sex)])
     legend(xy.coords(-3.5,7),pch= pch,col=colours,legend=smpls,bty="n",cex = 0.5) # rep(c(20,17,length(pc$x[,3])),8)
     mar=c(3.1,2.1,3.1,0.1)
     par(mar=mar)
   }
   
   
   out2<-combn(var.list,2)
   mapply(pcrld.plot,out2[1,],out2[2,],USE.NAMES=FALSE)
   
 }
 
 
 
 
 
 
 
 
 ### ==============================
 ## arrayQualityMetrics
 ### ==============================
 message("The array quality metrics is under procedure")
 suppressMessages(arrayQualityMetrics(ExpressionSet(assayData=vst),
                                      outdir=Metrics))
 
}

## Nico 20) Good job. It starting to have the right shape this function. I think it's good for you also to see the different steps. I.e. you start by collecting code, then you turn it into a function, then you wrap it into an executable script, then you optimise the function and finally integrate it into a package or another environment. That's good learning and shows you the different steps of a function conception from a prototype to a "production" stable function. As you'll develop more and more, you'll be able to skip some of these steps or at least streamline them, in the sense that you'll have gathered experience and knowledge and will simply design better prototypes. Keep on the good work, I enjoy it :-)

#genericRnaSeq()



Main <- function(){
  ### ================ main
  ## define the arguments
  option_list <- list(
    make_option(c("-op", "--output_prefix"),dest="op", type="character", default="",
                help="The output prefix, if wanted"),
    make_option(c("-s", "--samples"),dest="s", type="character", default="samples.csv",
                help="The filename for the samples file, if changed"),
    make_option(c("-wd", "--wd_template"), dest="wd", type="character", default="~/Git/UPSCb",
                help="The wd template, if other than default"),
    make_option(c("-hd", "--htseq_dir"),dest="hd", type="character", default="~/Git/UPSCb/HTSeq",
                help="The HTSeq directory, if other than default"),
    make_option(c("-hp", "--htseq_pattern"),dest="hp", type="character", default="*.txt",
                help="The HTSeq pattern, if other than default"),
    make_option(c("-hg", "--htseq_gsub"),dest="hg", type="character", default="_.*_STAR\\.txt",
                help="The HTSeq substitution (gsub), if other than default"),
    make_option(c("-hmc", "--htseq_mc_cores"),dest="hmc", type="integer", default=9,
                help="Number of mc cores for the HTSeq"),
    make_option(c("-r", "--reads"), dest="r", type="integer", default=2,
                help="Single or double reads (1 or 2)"),
    make_option(c("-p1", "--pca1"), dest="p1", type="integer", default=1,
                help="What dimension the first pca axis is going to show"),
    make_option(c("-p2", "--pca2"), dest="p2", type="integer", default=2,
                help="What dimension the second pca axis is going to show"),
    make_option(c("-p3", "--pca3"),dest="p3", type="integer", default=3,
                help="What dimension the third pca axis is going to show"),
    make_option(c("-m", "--metrics"),dest="m", type="character", default="Metrics",
                help="The output directory for the arrayQualityMetrics"),
    make_option(c("-pna", "--palette_name"),dest="pna", type="character", default="Dark2",
                help="The forward fastq file"),
    make_option(c("-pnr", "--palette_nr"),dest="pnr", type="integer", default=8,
                help="Number of different colors in the palette"),
    make_option(c("-nf", "--norm_fun"),dest="nf", type="character", default="vst",
                help="If not vst should be performed, please write rld"))
  
  ## get the options
  opt <- parse_args(OptionParser(option_list=option_list))
  
  genericRnaSeq(opt$op,opt$s,opt$wd,opt$hd,
                opt$hp,opt$hg,opt$hmc,opt$r,
                opt$p1,opt$p2,opt$p3,opt$m,
                brewer.pal(opt$pna,opt$pnr),opt$nf) 
  
  return("GenericRnaSeq completed")
  
}
Main()
#if(interactive()){
#  cat(Main()) Doesn't seem to work but otherwise useful
}
  

