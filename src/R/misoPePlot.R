## set the wd
setwd("/gulo/proj_nobackup/b2012243/data/accdata/MISO/PE_distribution")

## read the header
mat <- do.call(rbind,strsplit(sapply(dir(path=".",pattern="*.insert_len$",recursive=TRUE,full.names=TRUE),scan,what="character",n=1),"=|,"))

## mean
png(file="PE-fragment-size-distribution-boxplots.png",width=600,height=600,pointsize=16)
par(mfrow=c(2,2))
boxplot(as.numeric(mat[,2]),xlab="mean")

## sd
boxplot(as.numeric(mat[,4]),xlab="sd")

## dispersion
boxplot(as.numeric(mat[,6]),xlab="dispersion")

## number of reads used
boxplot(as.numeric(mat[,8]),xlab="number of reads")

dev.off()

