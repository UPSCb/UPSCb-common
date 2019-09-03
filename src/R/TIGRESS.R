suppressPackageStartupMessages({
    library(matrixStats)
    library(data.table)
    library(lars)
    library(parallel)
    library(Matrix)
    library(methods)
})

args <- commandArgs(trailingOnly = TRUE);
inf <- args[1]
outf <- args[2]
cores <- args[3]

tic <- Sys.time()

difft <- function(ti = tic){
      toc <- Sys.time()
      ti <- round(difftime(toc,ti, units="secs"))
      return(ti)
  }

statusUpd <- function(mess){
    message(paste("[\t",difft(),"]\t",mess))
}

TIGRESS <- function(in.file, out.file, cores=2){
  statusUpd("Reading file")
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
  statusUpd("Starting loop")
  mclapply(1:p,mc.cores = cores, function(j){
      tk <- Sys.time()
      s <- setdiff(seq(p),j)

      x <- mdat[,s]
      y <- mdat[,j,drop=FALSE]
      
      halfsize <- as.integer(n/2)

      freq <- sparseMatrix(nstepsLARS,p-1)
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
              r <- lars(xs[i1,],y[i1],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE,use.Gram=FALSE)
              r2 <- lars(xs[i2,],y[i2],max.steps=nstepsLARS,normalize=FALSE,trace=FALSE,use.Gram=FALSE)

                                        #test if lars produced valid outcome and record frequencies
              if(length(r$lambda)==nstepsLARS & length(r2$lambda)==nstepsLARS){
                  freq<-freq + abs(sign(r$beta[2:(nstepsLARS+1),]))
                  freq<-freq + abs(sign(r2$beta[2:(nstepsLARS+1),]))
                                        #only increment i if stability selection was successful in order to guarantee nBootstraps are actually
                                        #produced
                  i = i+1
              }
          }
      }

      freq <- freq/(2*nbootstrap)
      re <- colCumsums(as.matrix(freq))[5,]/5
      res <- data.table("i" = genes[j], "j" = genes[s], "x" = re)[x>0]
      tf <- tempfile(tmpdir = tdir)
      write.table(res,file = tf, quote = FALSE, row.names = FALSE,
                  col.names = FALSE, sep = "\t")
      statusUpd(paste("Finished gene",j,"time:",difft(tk)))
  })  
}

TIGRESS(inf, outf, cores)
