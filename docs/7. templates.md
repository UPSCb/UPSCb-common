# Templates

Before **re-inventing the wheel**, check the templates directory! A number of useful templates are available there:

1. This directory contains code templates for Biological quality analysis and Differetnial expression analysis in the RNA-Seq pipeline

    - BiologicalQA.R: a template to perform the biological QA of RNA-Seq data quantified by salmon
    - DifferentialExpression.R: a template to perform the differential expression of RNA-Seq data using DESeq2

    The templates contain `CHANGEME` tokens, where changes need to be made to adjust to your project. The `R` templates are written using the `#'` comments that allow for the knitting of `Rmd` files and `html` reports.

    Some of the `CHANGEME` tokens are within markdown code block. These are meant to be invisible when knitting (`eval=FALSE, echo=FALSE`), so you need *NOT* change them, just follow the instructions they contain to edit the subsequent code.

    The templates are not completely self-explanatory and as such the report generated from them will not be. Spend time to comment and interpret the results you observe using `#'` comments. These will be turned into text when knitting.

    The templates expect the following directory structure. 
    ```text
    project/  
    ├──  src/R        - where your custom R code resides
    ├──  doc          - where the relevant sample information can be found
    ├──  UPSCb-common - the checkout of the GitHub UPSCb/UPSCb-common repository, best as a submodule to your project
    ├── data         - a link to the data directory
    │    └── salmon - the directory containing the salmon results. One directory per sample, each containing the salmon results for that sample
    └── reference a link to the reference genomic information (_e.g._ transcript to gene mapping file)
    ```


2. Helper files

    These are styling files you need for proper formatting of the reports from these scripts as well as adding some logos etc.
    
    - bulogo2.png 
    - style.css: Controls the width and margins of the report
    - header.html: This adds the little Github logo in the top corner
    - footer.html: this adds the logo and the links at the bottom of the report. You need to edit for your name and links
    - Mfuzz.R: Only for co-expression analysis
