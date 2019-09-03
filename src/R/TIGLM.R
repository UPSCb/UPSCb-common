library(matrixStats)
library(data.table)
library(glmnet)
TIGLM <- function(in.file, out.file, cores=2){

  dat <- fread(in.file)
  mdat <- as.matrix(dat)
  out.d <- dirname(out.file)
  dir.create(paste(out.d,"/tmp",sep=""))
  tdir <- paste(out.d,"/tmp",sep="")
  sdat <- sign(mdat)
  sel <- colSums(sdat) > 5
  mdat <- mdat[,sel]
  nbootstrap <- 1000
  nstepsLARS <- 5
  alpha <- 0.2

  n <- nrow(mdat)
  p <- ncol(mdat)

  genes <- colnames(mdat)

  mclapply(1:p,mc.cores = cores, function(j){

    s <- setdiff(seq(p),j)

    x <- mdat[,s]
    y <- mdat[,j,drop=FALSE]

    halfsize <- as.integer(n/2)

    freq <- sparseMatrix(p-1,nstepsLARS)
    i <- 0

    while(i<nbootstrap) {
      # Randomly reweight each variable
      xs <- t(t(x)*runif(p-1,alpha,1))

      # Ramdomly split the sample in two sets
      perm <- sample(n)
      i1 <- perm[1:halfsize]
      i2 <- perm[(halfsize+1):n]

      #for lowly expressed samples it is possible that only 0 values will get selected
      #this would fail the lars
      if(! (sum(y[i1])==0 || sum(y[i2])==0) ){
        # run LARS on each randomized, sample and check which variables are selected
        r <- glmnet(xs[i1,],y[i1],nlambda=nstepsLARS,family="gaussian",
                    standardize = FALSE, dfmax = 50)
        r2 <- glmnet(xs[i2,],y[i2],nlambda=nstepsLARS,family="gaussian",
                     standardize = FALSE, dfmax = 50)

        #test if lars produced valid outcome and record frequencies
        if(length(r$lambda)==nstepsLARS & length(r2$lambda)==nstepsLARS){
          freq<-freq + abs(sign(r$beta))
          freq<-freq + abs(sign(r2$beta))
          #only increment i if stability selection was successful in order to guarantee nBootstraps are actually
          #produced
          i = i+1
        }
      }
    }

    freq <- freq/(2*nbootstrap)
    re <- rowSums(freq)/nstepsLARS
    res <- data.table("i" = genes[j], "j" = genes[s], "x" = re)
    tf <- tempfile(tmpdir = tdir)
    write.table(res,file = tf, quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = "\t")
  }
  )

}