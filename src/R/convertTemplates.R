#' ---
#' title: "Template conversion"
#' author: "Nicolas Delhomme"
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
  library(here)
  library(knitr)
})

#' # Rmd to R
#' 
#' ## Differential Expression
purl(here("template/R/DifferentialExpression_WithGOenrichment.Rmd"),
     output=here("template/R/DifferentialExpression_WithGOenrichment.R"),
     documentation=2)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
