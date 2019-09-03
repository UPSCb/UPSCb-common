suppressPackageStartupMessages(require(GenomicRanges))

args <- commandArgs(trailingOnly=TRUE);

inf <- args[1]
chr <- args[2]
st <- as.integer(args[3])
en <- as.integer(args[4])

load(inf)

gr <- grep(".GR$",ls(),value=TRUE)

target <- GRanges(chr,IRanges(st,en))
res <- sort(subsetByOverlaps(get(gr), target))

write.table(as.data.frame(res), sep="\t", quote=FALSE)
