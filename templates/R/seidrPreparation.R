#' ---
#' title: "CHANGEME Network data preparation"
#' author: "CHANGEME and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(here)
  library(readr)
})

#' Helper
source(here("UPSCb-common/src/R/featureSelection.R"))

#' # Data
#' ```{r CHANGEME2,eval=FALSE,echo=FALSE}
#'  CHANGEME is the variance stabilised data, where the transformation was done taking the
#'  model into account (_i.e._ `blind=FALSE`)
#' ```
vst <- read_tsv(here("CHANGEME")) %>% 
  column_to_rownames("ID")

#' # Filter
#' Different ways of filtering are considered. An cutoff on the absolute expression level
#' or on the variance (standard score), respectively. The latter might be more adequate to
#' avoid filtering lowly expressed genes that show sufficient, metadata relevant, variance.
sels <- rangeFeatureSelect(counts=as.matrix(vst),
                           conditions=dds$Tissue,
                           nrep=6)

e.cutoff <- 1

s.sels <- rangeFeatureSelect(counts=as.matrix(vst),
                             conditions=dds$Tissue,
                             nrep=6,scale=TRUE)
s.cutoff <- 1

table(s.sels[[s.cutoff]],sels[[e.cutoff  + 1]])
sum(featureSelect(as.matrix(vst),dds$Tissue,exp=0.1,6))
sum(featureSelect(as.matrix(vst),dds$Tissue,exp=0.05,6,scale=TRUE))

#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#'  CHANGEME is the vst cutoff devised from the plot above. The goal is to remove / reduce
#'  the signal to noise. Typically, this means trimming the data after the first sharp decrease 
#'  on the y axis, most visible in the non logarithmic version of the plot
#' ```
vst.cutoff <- "CHANGEME"

#' # Export
dir.create(here("data/seidr"),showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(vst[sels[[vst.cutoff+1]],]),
            file=here("data/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.[0-9]+$","",rownames(vst)[sels[[vst.cutoff+1]]])),
            file=here("data/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
