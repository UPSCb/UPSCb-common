"percentile" <- function(x,probs=seq(0,1,.01),...){
  quantile(x,probs=probs)
}