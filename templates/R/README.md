# Templates

This directory contains code templates:

1. RNA-Seq
  1. BiologicalQA.R: a template to perform the biological QA of RNA-Seq data quantified by salmon
  2. DifferentialExpression.R: a template to perform the differential expression of RNA-Seq data using DESeq2

The templates contain `CHANGEME` tokens, where changes need to be made to adjust to your project. The `R` templates are written using the `#'` comments that allow for the knitting of `Rmd` files and `html` reports.

Some of the `CHANGEME` tokens are within markdown code block. These are meant to be invisible when knitting (`eval=FALSE, echo=FALSE`), so you need *NOT* changing them, just follow the instructions they contain to edit the subsequent code.

The templates are not completely self-explanatory and as such the report generated from them will not be. Spend time to comment and interpret the results you observe using `#'` comments. These will be turned into text when knitting.

## RNA-Seq

The templates expect the following directory structure. 

  - project -
            |- src/R        - where your custom R code resides
            |- doc          - where the relevant sample information can be found
            |- UPSCb-common - the checkout of the GitHub UPSCb/UPSCb-common repository, best as a submodule to your project
            |- data         - a link to the data directory
                  |- salmon - the directory containing the salmon results. One directory per sample, each containing the salmon results for that sample
            |- reference a link to the reference genomic information (_e.g._ transcript to gene mapping file)
            
