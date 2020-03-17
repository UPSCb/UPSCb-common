#' ---
#' title: "Seidr ROC curves"
#' author: "CHANGEME"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(gplots)
  library(here)
  library(matrixStats)
  library(pander)
  library(pracma)
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
  
  aucs <- unlist(dat %>% group_by(ALGORITHM) %>% group_map(.,function(x,...){round(trapz(x$FP,x$TP),digits=3)}))
  names(aucs) <- unique(dat$ALGORITHM)
  
  p <- ggplot(dat,aes(x=FP,y=TP,col=ALGORITHM,group=ALGORITHM)) +
    geom_line() + 
    scale_x_continuous(name="1-specificity (FPR)") + 
    scale_y_continuous(name="sensitivity (TPR)") +
    ggtitle(label=paste(sub("_roc\\.tsv","",basename(f)), " ROC curve"))
  
  suppressMessages(suppressWarnings(plot(p)))
  
  return(list(stats=tabs,auc=aucs))
  
}

#' # ROC
#' ## Aggregated
#' ```{R CHANGEME1, echo=FALSE, eval=FALSE}
#' Change the path to the aggregated ROC results, if required. That ROC file must have been created using seidr roc -a option
#' ```
res <- plotRoc(here("data/seidr/roc/aggregated_roc.tsv"))

#' ### Stats of the gold standard analysis
pander(res$stats)

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
sts <- lapply(names(files),function(n,resb){
  resb[[n]]$stats
},resb)

names(sts) <- names(files)

pander(sts)

#' # Summary
#' Report all AUCs
aucs <- cbind(sapply(resb,"[[","auc"),aggregated=res$auc)

pander(aucs)

#' ## Heatmaps
#' ### AUCs
heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1))

heatmap.2(aucs,trace="none",col=hpal,margins=c(7.1,7.1),Colv=FALSE,dendrogram="row")

#' ### AUCs penalised by the number of Gold Standard (GS) edges used
#' The rationale here is to check the effect of a limited number of GS on the 
#' AUC calculation
resb$aggregated=res

gsNum <- sapply(lapply(lapply(resb,"[[","stats"),"[",2:3),rowSums)

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

