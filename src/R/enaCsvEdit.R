# Read a csv from the ENA submission and correct it
# a ; csv
dat <- read.csv2("~/Git/UPSCb/projects/spruce-needles-drought-stress/doc/ENA/ENA_Submission_PA_Drougth_JH_1.csv",
                 as.is=TRUE)

# order
dat <-dat[order(dat$SampleName),]

# look
str(dat)
length(unique(dat$SampleName))
length(unique(dat$SequencingDate))

# date
#dat$SequencingDate <- "2015-02-16T10:00:00"
table(dat$SequencingDate)
dat$SequencingDate <- ifelse(dat$SequencingDate=="2015_10_20","2015-10-20T10:00:00","2015-11-18T10:00:00")

# Sample Name
table(dat$SampleDescription)
table(dat$SampleName)
dat$SampleName <- sub("_[1,2]$","",dat$SampleName)
length(unique(dat$SampleName))
dat$SampleName <- paste(dat$SampleName,"biological replicate",rep(rep(1:3,each=4),6),"technical replicate",rep(rep(1:2,each=2),24))

# Description
length(unique(dat$SampleDescription))
length(unique(paste(dat$SampleDescription,"biological replicate",rep(rep(1:3,each=4),6),"technical replicate",rep(rep(1:2,each=2),24))))
dat$SampleDescription <- paste(dat$SampleDescription,"biological replicate",rep(rep(1:3,each=4),6),"technical replicate",rep(rep(1:2,each=2),24))

# save
write.csv(dat,file="~/Git/UPSCb/projects/spruce-needles-drought-stress/doc/ENA/ENA_Submission_PA_Drougth_JH_1.csv",
row.names=FALSE,quote=FALSE)
