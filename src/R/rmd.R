rmd <- function(x){
  mean(abs(x-mean(x)))/mean(x)
}
