library(minet)
setwd(args[2])
args <- commandArgs(trailingOnly = TRUE)
exp.dat <- fread(args[1])
mim <- build.mim(as.matrix(exp.dat))
ar <- aracne(mim)
arl <- ar[lower.tri(ar)]

gn <- colnames(ar)
fout <- "predictions.txt"
arli <- 1

for(i in 2:nrow(ar)){
  for(j in 1:(i-1)){
    if(arl[arli] > 0){
      cat(gn[i],gn[j],arl[arli],"\n",append=TRUE,file=fout,sep="\t")
    }
    arli <- arli + 1
  }
}