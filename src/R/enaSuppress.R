#' function
suppress <- function(id){
  message(paste(
    "<ACTION>",
    paste0('<SUPPRESS target="',id,'"/>'),
    "</ACTION>\n"))
}

echo.commands <- function(submission.filename){
  message( sprintf('curl -k -F "SUBMISSION=@%s" "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%%20Webin-763%%20AdDfb2VJ"',submission.filename))
  message( sprintf('curl -F "SUBMISSION=@%s" "https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%%20Webin-763%%20AdDfb2VJ"',submission.filename))
}

#' # UPSC-0022
#' ERX1659303 to ERX1659306
#' ERR1588657 to ERR1588660
#' ERS1294186 to ERS1294189
samples <- sprintf("ERS129418%d",6:9)
experiments <- sprintf("ERX165930%d",3:6)
runs <- sprintf("ERR15886%d",57:60)
suppress(c(samples,experiments,runs))
echo.commands("UPSC-0022-Suppress.Submission.xml")
