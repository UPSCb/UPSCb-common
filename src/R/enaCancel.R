#' function
cancel <- function(id){
  message(paste(
    "<ACTION>",
    paste0('<CANCEL target="',id,'"/>'),
    "</ACTION>\n"))
}

echo.commands <- function(submission.filename){
  message( sprintf('curl -k -F "SUBMISSION=@%s" "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%%20Webin-763%%20AdDfb2VJ"',submission.filename))
  message( sprintf('curl -F "SUBMISSION=@%s" "https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%%20Webin-763%%20AdDfb2VJ"',submission.filename))
}

#' # UPSC-0095
#' IDs
samples <- c("ERS1505219","ERS1507800","ERS1505213","ERS1505222","ERS1506707","ERS1507806","ERS1505226","ERS1507817","ERS1505233","ERS1507820","ERS1505236","ERS1507823","ERS1505293","ERS1508477")
cancel(samples)

experiments <- c("ERX1859050","ERX1864341","ERX1859044","ERX1859053","ERX1862394","ERX1864344","ERX1859056","ERX1864354","ERX1859063","ERX1864357","ERX1859066","ERX1864360","ERX1859144","ERX1865606")
cancel(experiments)

#' # UPSC-0128
#' ERS1801496 to ERS1801499; ERX2077399 to ERX2077402; ERR2017817 to ERR2017820
samples <- sprintf("ERS180149%d",6:9)
experiments <- sprintf("ERX2077%d",399:402)
runs <- sprintf("ERR20178%d",17:20)
cancel(c(samples,experiments,runs))
echo.commands("UPSC-0128.Submission.xml")
