#' ---
#' title: "Differential Expression"
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
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' TODO remember to add a function to check for correlation between logFC and 
#' transcript length - check for the effective length

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values accross your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=MDay,y=value,col=MGenotype,group=MGenotype)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter)
    } else {
        res <- results(dds,contrast=contrast,filter = filter)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
      message(sprintf(paste(
        ifelse(sum(sel)==1,
               "There is %s gene that is DE",
               "There are %s genes that are DE"),
        "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
        sum(sel),padj,
        lfc,expression_cutoff))
    }
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
        if (sum(val) >0){
          warning(sprintf(paste(
            ifelse(sum(val)==1,
                   "There is %s DE gene that has",
                   "There are %s DE genes that have"),
            "no vst expression in the selected samples"),sum(val)))
          sel[sel][val] <- FALSE
        } 

        if(export){
            if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
            write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot & sum(sel)>1){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=FALSE,
                                     plot=TRUE,export_pdf=TRUE,
                                     default_dir=here("analysis/DE"),
                                     default_prefix="DE"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))

    # sanity
    if(is.null(unlist(enrichment)) | length(unlist(enrichment)) == 0){
        message("No GO enrichment for",names(enrichment))
    } else {
        # write out
        if(export){
            write_tsv(bind_rows(enrichment),
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment$id,url=url,mc.cores=16L),
                    file=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
      if(plot){
        gocatname <- c(BP="Biological Process",
                       CC="Cellular Component",
                       MF="Molecular Function")
        degname <- c(all="all DEGs",
                     up="up-regulated DEGs",
                     dn="down-regulated DEGs")
        lapply(names(enrichment),function(n){
          lapply(names(enrichment[[n]]),function(de){
            lapply(names(enrichment[[n]][[de]]),function(gocat){
              dat <- enrichment[[n]][[de]][[gocat]]
              if(is.null(dat)){
                message("No GO enrichment for ",degname[de]," in category ",gocatname[gocat])
              } else {
                dat$GeneRatio <- dat$Significant/dat$Annotated
                dat$adjustedPvalue <- parse_double(sub("<","",dat$FDR))
                stopifnot(all(! is.na(dat$adjustedPvalue)))
                dat$Count <- as.numeric(dat$Significant)
                dat <- dat[order(dat$GeneRatio),]
                dat$Term <- factor(dat$Term, levels = unique(dat$Term))
                p <- ggplot(dat, aes(x =Term, y = GeneRatio, color = adjustedPvalue, size = Count)) + 
                  geom_point() +
                  scale_color_gradient(low = "red", high = "blue") +
                  theme_bw() + 
                  ylab("GeneRatio") + 
                  xlab("") + 
                  ggtitle(label=paste("GO enrichment:",degname[de]),
                          subtitle=(gocatname[gocat])) +
                  coord_flip()
                
                print(p)
                
                ggsave(file=here(default_dir,
                                 paste(default_prefix,n,de,gocat,"genes_GO-enrichment.pdf",sep="_")),
                       p)
              }
            })
          })
        })
      }
    }
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/salmon/dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' ## Gene of interests
#' ```{r goi, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you can plot the expression pattern of your gene of
#' interest. You need to have a list of genes in a text file, one geneID per line
#' The ID should exist in your vst data.
#' Note that for the plot to work, you also need to edit the first function (`line_plot`)
#' at the top of this file
#' ```
goi <- read_lines(here("doc/goi.txt"))
stopifnot(all(goi %in% rownames(vst)))
dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
#' #' ```{r res, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to define the contrast you want to 
#' study - see the example in the next block. 
#' 
#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' 
#' The `label` argument is typically one (or a combination) of the metadata stored in colData
#' 
#' The function allows for looking at the independent filtering results using `debug=TRUE`
#' 
#' If you are not satisfied with the default from DESeq2, you can set your own cutoff using `expression_cutoff`
#' 
#' You can change the default output file prefix using `default_prefix`
#' 
#' You can select the set of samples to be added to the `heatmap`, using the `sample_sel` argument. It takes a logical vector.
#' 
#' ```

#' ```{r contrast, echo=FALSE,eval=FALSE}
#'contrast1 <- extract_results(dds=dds,vst=vst,contrast="CHANGEME")
#' ```

#' ```{r sanitycheck, echo=FALSE,eval=FALSE}
#' CHANGEME - Here is the sanity check for complex contrasts in DE
#' And this is optional, you don't have to run this if you're not in doubt with the contrasts.
#' 
#' In the examples below, we assume that we wanted to compare between two subsets from the whole dataset
#' First chunk will be by using complex contrasts on the whole dataset. This is what we will use for the final result.
#' Second chunk is subsetting the dataset then apply simpler contrast. This will reduce the power (p-values will change).
#' We check whether the log2FC of the same gene are similar by plotting a scatter plot.
#' 1. If yes (y=x; some jitters allowed), then we did not make any mistake on setting complex contrast, we can continue to use the result from the first chunk.
#' 2. If no (cloud of points), we need to review if we did anything wrong on contrast specification in the first chunk.
#' 3. If yes (y=x) but there is a few points that is very far from y=x line, then it might come from accumulation of rounding error in control samples.
#' If we encounter rounding error accumulation, we should use the whole dds but relevel the ref so that the ref is within the pair that we want to compare.
#' 
#' Here assuming we have one control, two treatments, and 2 types of tissue in each condition. Tissue#1 is the control.
#' Then, for this specific contrast, we want to compare tissue#2 between two treatments conditions.
#' ```

#' ```{r complexcontrast, echo=FALSE,eval=FALSE}
#'resultsNames(dds)
#'tissue2_tr1tr2 <- extract_results(dds=dds,vst=vst,contrast=c(0,0,1,-1,1,-1),
#'.                                 default_dir = here("data/analysis/DE/test"),
#'.                                 default_prefix = "tissue2_tr1tr2_contrast_",
#'.                                 labels=paste(dds$Treatment,dds$Tissue,sep="_"),
#'.                                 sample_sel = dds$Tissue == "2" & dds$Treatment %in% c("Tr1","Tr2"), 
#'.                                 cexCol=.9, plot = TRUE, verbose = TRUE)
#' ```

#' ```{r simplercontrast, echo=FALSE,eval=FALSE}
#'sample_sel <- dds$Tissue == "2" & dds$Treatment %in% c("Tr1","Tr2")
#'ddssubset <- dds[,sample_sel]
#'ddssubset$Treatment <- relevel(ddssubset$Treatment,ref = "Tr1")
#'design(ddssubset) <- ~ Treatment
#'ddssubset$Treatment <- droplevels(ddssubset$Treatment)
#'ddssubset$Tissue <- droplevels(ddssubset$Tissue)
#'ddssubset <- DESeq(ddssubset)
#'vsdsubset <- varianceStabilizingTransformation(ddssubset,blind=FALSE)
#'vstsubset <- assay(vsdsubset)
#'vstsubset <- vstsubset - min(vstsubset)
#'resultsNames(ddssubset)
#'tissue2_tr1tr2 <- extract_results(dds=ddssubset,vst=vstsubset,contrast="Treatment_Tr2_vs_Tr1",
#'                                  default_dir = here("data/analysis/DE/test"),
#'                                  default_prefix = "tissue2_tr1tr2_subset_",
#'                                  labels=paste(ddssubset$Treatment,ddssubset$Tissue,sep="_"),
#'                                  cexCol=.9, plot = TRUE, verbose = TRUE)
#'contrast_result <- na.omit(read.csv(here("data/analysis/DE/test/tissue2_tr1tr2_contrast_results.csv")))
#'subset_result <- na.omit(read.csv(here("data/analysis/DE/test/tissue2_tr1tr2_subset_results.csv")))
#'colnames(contrast_result) <- paste0("con_",colnames(contrast_result))
#'contrast_result$gene <- contrast_result$con_X
#'colnames(subset_result) <- paste0("sub_",colnames(subset_result))
#'subset_result$gene <- subset_result$sub_X
#'dat <- inner_join(contrast_result,subset_result)
#'ggplot(dat, aes(x=con_log2FoldChange, y=sub_log2FoldChange)) + 
#'  geom_point() + 
#'  labs(title = "Tissue#2 Tr2 vs Tr1")
#'remove(ddssubset,vstsubset,vsdsubset)
#' ```

#' ### Venn Diagram
#' ```{r venn, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you typically would have run several contrasts and you want
#' to assess their overlap plotting VennDiagrams.
#' 
#' In the examples below, we assume that these resutls have been saved in a list
#' called `res.list`
#' ```

#' #### All DE genes
#' ```{r venn2, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","all"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' #### DE genes (up in mutant)
#' ```{r venn3, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","up"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' #### DE genes (up in control)
#' ```{r venn4, echo=FALSE,eval=FALSE}
#'grid.newpage()
#'grid.draw(venn.diagram(lapply(res.list,"[[","dn"),
#'                       NULL,
#'                       fill=pal[1:3]))
#' ```

#' ### Gene Ontology enrichment
#' ```{r go, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table.
#' We used to export another table consisting
#' of only the `id` and `padj` columns for using as input for _e.g._
#' REVIGO; but since flash is EOL and REVIGO not updated, we instead rely on 
#' the RtoolBox treemap.
#' 
#' In addition we now also export the list of genes that most likely resulted in
#' the corresponding go enrichment.
#' 
#' Make sure to change the `url` to match your species
#' 
#' ```
#' TODO USE THE independent filtering to decide on the background. Think about it
background <- rownames(vst)[featureSelect(vst,dds$MGenotype,exp=CHANGEME)]

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="athaliana")
})

dev.null <- lapply(names(enr.list),function(n){
    lapply(names(enr.list[[n]]),function(de){
        extractEnrichmentResults(enr.list[[n]][[de]],
                                 diff.exp=de,
                                 genes=res.list[[n]][[de]],
                                 default_prefix=paste(n,de,sep="-"),
                                 url="athaliana")
    })
})

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


