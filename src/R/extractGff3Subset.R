#' ---
#' title: "Extract Gff3 subtract"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("~/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/")
#' ```

#' Libraries
suppressPackageStartupMessages(library(genomeIntervals))

#' Source the helper file
source("~/Git/UPSCb/src/R/gff3Utilities.R")

#' Read the gff3 file
gff3 <- readGff3("/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3",
                 quiet=TRUE)

#' Define the gene list
gene.list <- c("MA...","MA...")

#' # Extract
#' The subset of interest
subgff3 <- extractFromGff3UsingGeneIDs(gff3=gff3,IDs=gene.list)

#' And save it
writeGff3(subgff3,file="A..FILE..NAME")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
#' 
