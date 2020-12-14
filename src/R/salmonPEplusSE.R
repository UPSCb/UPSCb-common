#' ---
#' title: "Salmon SE and PE data"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  #  library(data.table)
  #  library(DESeq2)
  #  library(gplots)
  library(here)
  #  library(hyperSpec)
  #  library(parallel)
  #  library(pander)
  #  library(plotly)
  library(LSD)
  library(tidyverse)
  library(tximport)
  #  library(vsn)
})

#' tx2gene translation table
tx2gene <- suppressMessages(read_delim("/mnt/picea/storage/reference/Populus-tremula/v2.2/annotation/tx2gene.tsv",delim="\t",
                                       col_names=c("TXID","GENE")))

#' # Raw data
filelist <- list.files(here("../../test"), 
                       recursive = TRUE, 
                       pattern = "quant.sf",
                       full.names = TRUE)

#' name the file list vector
names(filelist) <- basename(dirname(filelist))

#' Read the expression at the gene level
counts <- suppressMessages(round(tximport(files = filelist, 
                                          type = "salmon",
                                          tx2gene=tx2gene)$counts))

#' How many reads?
colSums(counts)

#' How do they correlate?
cor.test(counts[,1],counts[,2])

#' Plot them
#' 
#' Definitely the unpaired data are representative of the PE data
#' We should keep that in mind for dataset where the sequencing depth
#' might be limiting and simply consider the unpaired as a tech rep and
#' sum them up with the paired data to increase coverage
heatscatter(counts[,"PE"],counts[,"SE"],log="xy")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
