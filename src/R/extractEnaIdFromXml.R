#' # Libraries
library(tools)
#' # Functions
#' ## Compare the XML and CSV
"extractIDs" <- function(sample.xml=character(0),asmb.name){
    stopifnot(file.exists(sample.xml))
    
    # get the sample info and study ID
    info <- data.frame(
        SNAM=system(paste("grep TITLE",
                              sample.xml,
                              '| sed "s:\\<:\\>:g" | awk -F\\>',
                              "'{print $3}'"),
                        intern=TRUE),
        DESC=system(paste("grep DESCRIPTION",
                               sample.xml,
                               '| sed "s:\\<:\\>:g" | awk -F\\>',
                               "'{print $3}'"),
                         intern=TRUE),
        SaREF=system(paste("grep alias",
                             sample.xml, '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE)
    )
    
    # extract the study id
    info$StRef=sub("-S000[1-3]","",info$SaREF)
    info$fFile=asmb.name
    info$cFile=sub("fasta","list",info$fFile)
    info$fMd5 <- md5sum(info$fFile)
    info$cMd5 <- md5sum(info$cFile)
    
    # return
    return(info)
}

#' process
setwd("~/Git/UPSCb/projects/spruce-diversity/doc")

#' ## UPSC-0025, 59
i<-25
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- extractIDs(sample.xml,
           c("Pabies19540267.fasta.gz","Pabies19830406.fasta.gz"))

#' ## UPSC-0026
i<-26
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pwilsonii1963848.fasta.gz")))

#' ## UPSC-0027, 62
i<-27
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PasperataLiujq-09xz-lzt-190.fasta.gz",
            "Pasperata19800139.fasta.gz")))

#' ## UPSC-0028, 66
i<-28
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,
            c("Pcrassifolia827.fasta.gz","Pcrassifolia818.fasta.gz")))

#' ## UPSC-0029, 67
i<-29
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,
                          c("Pengelmannii5915.fasta.gz",
                            "Pengelmannii19780230.fasta.gz",
                            "Pengelmannii19780241.fasta.gz")))

#' ## UPSC-0030, 69
i<-30
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pglauca10755.fasta.gz",
            "Pglauca19490210.fasta.gz",
            "Pglauca19820049.fasta.gz")))

#' ## UPSC-0031, 70
i<-31
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,
                          c("Pglehnii13988.fasta.gz",
                            "Pglehnii19760375.fasta.gz")))

#' ## UPSC-0032,51,71
i<-32
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,
            c("Pjezoensis14410.fasta.gz",
              "Pjezoensis19760388.fasta.gz")))

i<-51
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pjezoensis14247.fasta.gz")))

#' ## UPSC-0033, 72
i<-33
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PkoraiensisZR09040.fasta.gz","Pkoraiensis19921086.fasta.gz")))

#' ## UPSC-0034, 74
i<-34
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Plikiangensis2684.fasta.gz",
  "Plikiangensis1-574.fasta.gz",
  "Plikiangensis19860201.fasta.gz")))

#' ## UPSC-0035, 75
i<-35
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pmariana19680218.fasta.gz","Pmariana19930565.fasta.gz")))

#' ## UPSC-0036, 77
i<-36
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pmexicana1991052.fasta.gz","Pmexicana19880558.fasta.gz")))

#' ## UPSC-0037, 78
i<-37
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PmeyeriO60923.fasta.gz","Pmeyeri19860107.fasta.gz")))

#' ## UPSC-0038, 79
i<-38
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pneoveitchiilzt-111.fasta.gz","PneoveitchiiTaiBai.fasta.gz")))

#' ## UPSC-0039, 81
i<-39
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Porientalis19830568.fasta.gz","Porientalis19870081.fasta.gz")))

#' ## UPSC-0040, 82
i<-40
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Ppungens19780256.fasta.gz","Ppungens19780259.fasta.gz")))

#' ## UPSC-0041, 83
i<-41
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pretroflexa2004.fasta.gz","Pretroflexa1872.fasta.gz")))

#' ## UPSC-0042, 84
i<-42
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Prubens.fasta.gz","Prubens6169.fasta.gz")))

#' ## UPSC-0043, 85
i<-85
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pschrenkiana15702.fasta.gz","Pschrenkiana1990060.fasta.gz",
  "Pschrenkiana19860197.fasta.gz")))

#' ## UPSC-0044, 86
i<-44
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Psmithiana19840137.fasta.gz","Psmithiana19860408.fasta.gz")))

#' ## UPSC-0045, 61
i<-45
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Palcoquiana19770679.fasta.gz")))

#' ## UPSC-0046, 60
i<-46
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PaurantiacaF012.fasta.gz")))

#' ## UPSC-0047, 63
i<-47
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pbrachytyla19930339.fasta.gz")))

#' ## UPSC-0048, 64
i<-48
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pbreweriana19760302.fasta.gz")))

#' ## UPSC-0049, 65
i<-49
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pchihuahuana19880553.fasta.gz")))

#' ## UPSC-0050, 68
i<-50
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pfarreri.fasta.gz")))

#' ## UPSC-0052, 73
i<-52
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pkoyamae24-670.fasta.gz")))

#' ## UPSC-0053, 76
i<-53
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pmaximowiczii18064.fasta.gz")))

#' ## UPSC-0054, 93
i<-54
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pobovata.fasta.gz")))

#' ## UPSC-0055
i<-55
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PpurpureaLiuJQ-QTP-2011-294.fasta.gz")))

#' ## UPSC-0056, 87
i<-56
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("PspinulosaF022.fasta.gz")))

#' ## UPSC-0057
i<-57
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Ptorano99251086.fasta.gz")))

#' ## UPSC-0058, 89
i<-58
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pwilsonii19800115.fasta.gz")))

#' ## UPSC-0080
i<-80
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pchihuahuana18063.fasta.gz")))

#' ## UPSC-0090
i<-90
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pomorika37358.fasta.gz")))

#' ## UPSC-0091
i<-91
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Psitchensis5294.fasta.gz")))

#' ## UPSC-0092
i<-92
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
df <- rbind(df,extractIDs(sample.xml,c("Pbrachytyla18159.fasta.gz")))

stopifnot(all(!duplicated(df$fFile)))
stopifnot(all(dir(".",pattern="*.fasta.gz") %in% df$fFile))

save(df,file="assembly.rda")
