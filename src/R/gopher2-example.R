# dat = your list of genes of interest (a subset of a population)
dat <- scan("~/Git/UPSCb/projects/facility/doc/tiggy-gene-example.txt",what="character")

# bg = your population (think what defines it.)
bg <- scan("~/Git/UPSCb/projects/facility/doc/Tiggy.spruce.background.txt",what="character",skip=1)

# we silence warnings (not good practice, but we know what we're doing) - think to adjust the path
suppressPackageStartupMessages(source("~delhomme/Git/UPSCb-common/src/R/gopher.R"))

# we just quantify the run time
# task has to be a list
# task can take any value from: go, kegg, mapman and pfam
system.time(enrichment <- gopher(dat,task = list("go","kegg","pfam"),background = bg,url="pabies"))

# if no background (i.e. the whole population is to be used)
system.time(enrichment <- gopher(dat,task = list("go","kegg","pfam"),url="pabies"))

# enrichment will contain a list of tibbles (a "type of" data.frame), on per "task"

# for go, you can for example export the GO ID and the FDR to a file and then 
# upload that to REVIGO (http://revigo.irb.hr)