suppressPackageStartupMessages(library(optparse))

Main <- function(){
  ### ================ main
  ## define the arguments
  option_list <- list(
    make_option(c("-op", "--output_prefix"),dest="op", type="character", default="",
                help="The output prefix, if wanted"))
  opt <- parse_args(OptionParser(option_list=option_list))
  
  return(opt$op)
}  
Main()
