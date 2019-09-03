args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
library(GeneNet)
library(data.table)
dat <- fread(args[2])
out <- args[3]

pcor.dyn = ggm.estimate.pcor(as.matrix(dat), method = "dynamic")
write.table(abs(pcor.dyn), out, quote = FALSE,
            col.names = FALSE, row.names = FALSE, sep = "\t")
