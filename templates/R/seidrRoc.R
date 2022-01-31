#' ---
#' title: "Seidr ROC curves"
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
  library(gplots)
  library(here)
  library(matrixStats)
  library(pander)
  library(tidyverse)
})

#' * Colors
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Functions
plotRoc <- function(f){
  
  dat <- read_tsv(f,
           col_names = c("TP","FP","PR","ALGORITHM"),
           col_types = cols(.default=col_double(),
                            ALGORITHM = col_character()),
           comment="#")
  
  vals <- read_lines(f) %>% subset(grepl("#",.))
  lvals <- length(vals)/3
  tabs <- tibble(ALGO=sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",2),
         TP=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",2)),
         FP=as.integer(sapply(strsplit(vals[1:lvals],"\t"),"[",3)),
         AUC=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(1,length.out=lvals,by=2)],"\t"),"[",1))),
         AUPR=as.double(sub(".* ","",sapply(strsplit(vals[lvals + seq(2,length.out=lvals,by=2)],"\t"),"[",1))),
         )
  
  # we used to calc the AUC before seidr did it, so that's left in as 
  # an example how to use trapz
  # require(pracma)
  # aucs <- dat %>% group_by(ALGORITHM) %>% 
  #   group_map(.,function(x,...){data.frame(ALGORITHM=unique(x$ALGORITHM),
  #                                          AUC=round(trapz(x$FP,x$TP),digits=3))},
  #             .keep=TRUE) %>% bind_rows()
  # 
  p <- ggplot(dat,aes(x=FP,y=TP,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="1-specificity (FPR)") + 
    scale_y_continuous(name="sensitivity (TPR)") +
    ggtitle(label=paste(sub("_roc\\.tsv","",basename(f)), " ROC curve"))
  
  suppressMessages(suppressWarnings(plot(p)))
  
  return(tabs)
  
}

#' # ROC
#' ## Aggregated
#' ```{R CHANGEME1, echo=FALSE, eval=FALSE}
#' Change the path to the aggregated ROC results, if required. That ROC file must have been created using seidr roc -a option
#' ```
res <- plotRoc(here("data/seidr/roc/aggregated_roc.tsv"))

#' ### Stats of the gold standard analysis
pander(res)

#' ## Backbone
#' ```{R CHANGEME1, echo=FALSE, eval=FALSE}
#' Change the path to the backbone ROC results directory as well as the file matching patter,if required. 
#' These ROC files must have been created using seidr roc -a option
#' ```
files <- dir(here("data/seidr/roc"),pattern="backbone.*\\.tsv",full.names=TRUE)
names(files) <- gsub("backbone-|_roc.tsv","",basename(files))
files <- files[order(as.integer(sub("-percent","",names(files))))]
resb <- lapply(files,plotRoc)

#' ### Stats of the gold standard (GS) analysis
pander(resb)

#' # Summary
#' Report all AUCs
aucs <- c(lapply(resb,select,c("ALGO","AUC")),aggregated=list(select(res,c("ALGO","AUC")))) %>% 
  enframe %>% unnest(cols=c("value")) %>% pivot_wider(names_from=name,values_from=AUC) %>% 
  column_to_rownames("ALGO") %>% as.matrix()

pander(aucs)

#' ## Heatmaps
#' ### AUCs
#' **NOTE** In the heatmap below, the only fair comparisons are to look
#' across the `aggregated` column or the `irp` row. Any other comparisons,
#' albeit possibly informative will be biased by the selection of a subset 
#' of edges for a given algorithm, making them theoretically incomparable.
#'
heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1))

heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' ### AUCs penalised by the number of Gold Standard (GS) edges used
#' The rationale here is to check the effect of a limited number of GS on the 
#' AUC calculation
resb$aggregated=res

gsNum <- lapply(resb,select,c("ALGO","TP","FP")) %>%
  enframe %>% unnest(cols=c("value")) %>% 
  mutate(TOTAL=TP+FP) %>% select(-c("TP","FP")) %>% 
  pivot_wider(names_from=name,values_from=TOTAL) %>% 
  column_to_rownames("ALGO") %>% as.matrix()

paucs <- aucs * t(t(gsNum) / colMaxs(gsNum))

pander(paucs)

heatmap.2(paucs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' # Conclusion
#' CHANGEME for some conclusion
#' ```{r empty, eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

