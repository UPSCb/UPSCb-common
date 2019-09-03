mailR <- function(to,subject,msg) {
  library(stringr)
  if(any(sapply(c("to","subject","msg"),function(f){class(get(f))}) != "character")) {
    stop("Please supply only character arguments to this function")
  }
  reg.test <- str_detect(to,"[A-Za-z0-9\\.]+\\@[A-Za-z0-9]+\\.[A-Za-z0-9\\.]+")
  if(!reg.test){
    stop("Your email doesn't appear to be valid")
  }
  if(nchar(subject)>2000){
    stop("Please restrict your subject to 2000 or less characters")
  }
  if(as.integer(object.size(msg))/1024^2>20){
    stop("Please restrict your message to be 20 MB or less")
  }
  to <- gsub('"',"'",to)
  subject <- gsub('"',"'",subject)
  msg <- gsub('"',"'",msg)
  cmd <- paste('echo','"',msg,'"','| mail -s "',subject,'" -a "From: RStudio"','--to "',to,'"')
  system(cmd)
}