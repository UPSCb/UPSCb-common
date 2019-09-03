#' # Functions
#' ## Compare the XML and CSV
"compareXmlToCsv" <- function(experiment.xml=character(0),
                              run.xml=character(0),
                              sample.xml=character(0),
                              sample.csv=character(0),
                              nfile=2){
    stopifnot(file.exists(experiment.xml))
    stopifnot(file.exists(run.xml))
    stopifnot(file.exists(sample.xml))
    stopifnot(file.exists(sample.csv))

    # get the experiment fields
    lst <- list(ES=data.frame(
        EID=rep(system(paste("grep alias",
                             experiment.xml,
                             '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE),
                each=nfile),
        ESID=rep(system(paste("grep SAMPLE_DESCRIPTOR",
                             experiment.xml,
                             '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE),
                each=nfile),
        LIB=rep(system(paste("grep LIBRARY_NAME",
                             experiment.xml,
                             '| sed "s:<:>:g" | awk -F\\>',
                             "'{print $3}'"),
                       intern=TRUE),
                each=nfile),
        SNAM=rep(system(paste("grep TITLE",
                              sample.xml,
                              '| sed "s:<:>:g" | awk -F\\>',
                              "'{print $3}'"),
                        intern=TRUE),
                 each=nfile),
        SDESC=rep(system(paste("grep DESCRIPTION",
                               sample.xml,
                               '| sed "s:<:>:g" | awk -F\\>',
                               "'{print $3}'"),
                         intern=TRUE),
                  each=nfile),
        SID=rep(system(paste("grep alias",
                             sample.xml, '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE),
                each=nfile)
    ),
    R = data.frame(RID=rep(system(paste("grep alias",
                                        run.xml, '| awk -F\\"', "'{print $2}'"),
                                  intern=TRUE),
                           each=nfile),
                   REID=rep(system(paste("grep refname",
                                         run.xml, '| awk -F\\"', "'{print $2}'"),
                                   intern=TRUE),
                            each=nfile),
                   FNAM=system(paste("grep filename",
                                     run.xml, '| awk -F\\"',
                                     "'{print $2}'",
                                     '| sed "s:.*/::"'),
                               intern=TRUE))
    )

    # validate
    #stopifnot(all(lst$ES$ESID==lst$ES$SID))
    stopifnot(all(lst$ES$LIB==lst$ES$SNAM))

    # next validation
    stopifnot(all(lst$R$REID %in% lst$ES$EID))

    # no empty entries
    stopifnot(all(lst$ES$LIB != ""))
    
    # build the csv
    xml <- cbind(lst$ES[match(lst$R$REID,lst$ES$EID),],lst$R)

    # read the csv
    csv <- read.csv(sample.csv)

    # validate
    stopifnot(nrow(csv)==nrow(xml))

    # check sample names
    stopifnot(all(xml$SNAM %in% csv$SampleName))
    
    # check that sample name and run corresponds
    res <- NULL
    if(!identical(split(csv$FileName,csv$SampleName),
                        split(xml$FNAM,xml$SNAM))){
        xml <- xml[match(csv$FileName,xml$FNAM),]
        res <- cbind(csv,xml)
        res <- res[which(as.character(csv$SampleName) != as.character(xml$SNAM)),]
        res <- res[order(res$SampleName,res$FileName),
                   c("SampleName","SNAM","FileName","SID","EID","RID","FNAM")]
    }

    # return
    return(res)
}

#' ## Fix the CSV
"fixXML" <- function(old.experiment.xml=character(0),
                     old.run.xml=character(0),
                     #old.sample.xml=character(0),
                     sample.csv=character(0)){

    stopifnot(file.exists(old.experiment.xml))
    stopifnot(file.exists(old.run.xml))
    stopifnot(file.exists(sample.csv))

    # get the experiment fields
    lst <- list(E=data.frame(
        EID=system(paste("grep alias",
                             old.experiment.xml,
                             '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE),
        ESID=system(paste("grep SAMPLE_DESCRIPTOR",
                             old.experiment.xml,
                             '| awk -F\\"', "'{print $2}'"),
                       intern=TRUE),
        ESNAM=system(paste("grep LIBRARY_NAME",
                              old.experiment.xml,
                              '| sed "s:\\<:\\>:g" | awk -F\\>',
                              "'{print $3}'"),
                        intern=TRUE)),
        # S=data.frame(
        #     SID=system(paste("grep alias",
        #                      old.sample.xml, '| awk -F\\"', "'{print $2}'"),
        #                intern=TRUE),
        #     SNAM=system(paste("grep TITLE",
        #                       old.sample.xml,
        #                       '| sed "s:\\<:\\>:g" | awk -F\\>',
        #                       "'{print $3}'"),
        #                 intern=TRUE)),
        R = data.frame(
            RID=rep(system(paste("grep alias",
                                 old.run.xml, '| awk -F\\"', "'{print $2}'"),
                           intern=TRUE),
                    each=2),
            REID=rep(system(paste("grep refname",
                                  old.run.xml, '| awk -F\\"', "'{print $2}'"),
                            intern=TRUE),
                     each=2),
            FNAM=system(paste("grep filename",
                              old.run.xml, '| awk -F\\"',
                              "'{print $2}'",
                              '| sed "s:.*/::"'),
                        intern=TRUE))
    )

    # collect the info we need
    info <- cbind(
#        lst$S,
        lst$E,
        lst$R[match(unique(lst$E$EID),lst$R$REID),])

    # read the csv
    csv <- read.csv(sample.csv)

    # read the experiment
    expt <- scan(old.experiment.xml,what="character",sep="\n")

    # correct
    info$CSID <- info$ESID[match(csv[match(info$FNAM,csv$FileName),"SampleName"],
                                gsub("&quot;","",info$ESNAM))]

    info$CEID <- info[match(info$CSID,info$ESID),"EID"]

    # and sort
    info <- info[order(info$CEID),]

    # The experiment xml
    sel <- grep("<EXPERIMENT alias=",expt)
    expt[sel] <- sapply(1:nrow(info),function(i,f,s,e){
        sub(as.character(f$CEID[i]),
            as.character(f$REID[i]),e[s[i]])
    },info,sel,expt)

    # sel <- grep("SAMPLE_DESCRIPTOR",expt)
    # expt[sel] <- sapply(1:nrow(info),function(i,f,s,e){
    #     sub(as.character(f$ESID[i]),
    #         as.character(f$CSID[i]),e[s[i]])
    # },info,sel,expt)

    # # The sample xml
    # smpl <- scan(old.sample.xml,what="character",sep="\n")
    # sel <- grep("<SAMPLE alias=",smpl)
    # smpl[sel] <- sapply(1:nrow(info),function(i,f,s,e){
    #     sub(as.character(f$SID[i]),
    #         as.character(f$CSID[i]),e[s[i]])
    # },info,sel,smpl)

    # write out
    write(expt,file=sub(".xml$","-fix.xml",old.experiment.xml))
#    write(smpl,file=sub(".xml$","-fix.xml",old.sample.xml))

    # return the info
    return(info)
}

# If nfile=1 corrects for technical replicates
"correctRunExperimentIDs" <- function(experiment.xml,
                          run.xml,
                          sample.csv,
                          nfile=2){
  
  stopifnot(file.exists(experiment.xml))
  stopifnot(file.exists(run.xml))
  stopifnot(file.exists(sample.csv))
  
  lst <- list(ES=data.frame(
    EID=rep(system(paste("grep alias",
                         experiment.xml,
                         '| awk -F\\"', "'{print $2}'"),
                   intern=TRUE),
            each=nfile),
    ESID=rep(system(paste("grep SAMPLE_DESCRIPTOR",
                          experiment.xml,
                          '| awk -F\\"', "'{print $2}'"),
                    intern=TRUE),
             each=nfile),
    LIB=rep(system(paste("grep LIBRARY_NAME",
                         experiment.xml,
                         '| sed "s:<:>:g" | awk -F\\>',
                         "'{print $3}'"),
                   intern=TRUE),
            each=nfile)
  ),
  R = data.frame(RID=rep(system(paste("grep alias",
                                      run.xml, '| awk -F\\"', "'{print $2}'"),
                                intern=TRUE),
                         each=nfile),
                 REID=rep(system(paste("grep refname",
                                       run.xml, '| awk -F\\"', "'{print $2}'"),
                                 intern=TRUE),
                          each=nfile),
                 FNAM=system(paste("grep filename",
                                   run.xml, '| awk -F\\"',
                                   "'{print $2}'",
                                   '| sed "s:.*/::"'),
                             intern=TRUE))
  )
  
  # read the csv
  csv <- read.csv(sample.csv)
  
  # map the exp with the run
  lst$R$SampleName <- csv[match(lst$R$FNAM,csv$FileName),"SampleName"]
  lst$R <- cbind(lst$R,lst$ES[match(lst$R$SampleName,lst$ES$LIB),])
  
  # read the original run file
  run <- scan(run.xml,what="character",sep="\n")
  
  # select the lines to edit
  sel <- grep("<EXPERIMENT_REF refname=",run)
  
  if(nfile==2){
    # identify the fwd/rev
    dups <- duplicated(lst$R[,-grep("FNAM",colnames(lst$R))])
    
    # check that they are pairs
    stopifnot(sum(duplicated(lst$R[,-grep("FNAM",colnames(lst$R))])) == nrow(lst$R)/2)
    
    # check that they are ordered by pairs
    stopifnot(all(
      sub("_1\\.f.*q\\.gz$","",lst$R$FNAM[seq(1,nrow(lst$R),2)])
      ==sub("_2\\.f.*q\\.gz$","",lst$R$FNAM[seq(1,nrow(lst$R),2)+1])))
    
    # edit
    run[sel] <- sapply(1:nrow(lst$R[!dups,]),function(i){
      sub(lst$R[!dups,"REID"][i],lst$R[!dups,"EID"][i],run[sel[i]])
    })
  } else {
    # simply replace the experiment ID
    run[sel] <- sapply(1:nrow(lst$R),function(i){
      sub(lst$R[i,"REID"],lst$R[i,"EID"],run[sel][i])
    })
  }
  # write out
  write(run,file=sub(".xml$","-fix.xml",run.xml))
  
  # report
  invisible(lst)
}


#' # Validation
#' ## UPSC-0014
#' Results are OK
setwd("~/Git/UPSCb/projects/spruce-needle-diurnal-series/doc/ENA")
experiment.xml<-"UPSC-0014/UPSC-0014.Experiment.xml"
run.xml <- "UPSC-0014/UPSC-0014.Run.xml"
sample.xml <- "UPSC-0014/UPSC-0014.Sample.xml"
sample.csv <- "../Spruce-diurnals-time-series-samples_annotated_141205.csv"
U14 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0015
#' No swap, just a changed sample name
setwd("~/Git/UPSCb/projects/spruce-wood-time-series/doc/ENA/")
experiment.xml<-"UPSC-0015.Experiment.xml"
run.xml <- "UPSC-0015.Run.xml"
sample.xml <- "UPSC-0015.Sample.xml"
sample.csv <- "Spruce_time_serie_ENA.csv"
U15 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

# UPSC-0016 - verify
setwd("~/Git/UPSCb/projects/spruce-LMD/doc/ENA/")
experiment.xml<-"UPSC-0016.Experiment.xml"
run.xml <- "UPSC-0016.Run.xml"
sample.xml <- "UPSC-0016.Sample.xml"
sample.csv <- "Lignifying-xylem-of-Norway-spruce_Samples.csv"
U16 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0017
#' OK
setwd("~/Git/UPSCb/projects/ethylene-insensitive/doc/ENA")
experiment.xml<-"UPSC-0017.Experiment.xml"
run.xml <- "UPSC-0017.Run.xml"
sample.xml <- "UPSC-0017.Sample.xml"
sample.csv <- "Sample-file-ACCproject-ENA.csv"
U17 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0018
#' OK
setwd("~/Git/UPSCb/projects/spruce-SE-germinants/doc")
experiment.xml<-"UPSC-0018.Experiment.xml"
run.xml <- "UPSC-0018.Run.xml"
sample.xml <- "UPSC-0018.Sample.xml"
sample.csv <- "ENA-SE-germinants.csv"
compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0019
#' All samples swapped.
setwd("~/Git/UPSCb/projects/spruce-callus-culture/doc/ENA/")
experiment.xml<-"UPSC-0019.Experiment.xml"
run.xml <- "UPSC-0019.Run.xml"
sample.xml <- "UPSC-0019.Sample.xml"
sample.csv <- "ENA.csv"
U19 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' Fix it
res <- fixXML(experiment.xml,run.xml,sample.csv)
fix.xml<-paste0("UPSC-0019.Experiment-fix.xml")
U19 <- compareXmlToCsv(fix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0020
#' OK
setwd("~/Git/UPSCb/projects/swasp-rnaseq/doc/ENA")
experiment.xml<-"UPSC-0020.Experiment.xml"
run.xml <- "UPSC-0020.Run.xml"
sample.xml <- "UPSC-0020.Sample.xml"
sample.csv <- "SwAsp.csv"
compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0021
#' OK
setwd("~/Git/UPSCb/projects/wood/doc")
experiment.xml<-"UPSC-0021.Experiment.xml"
run.xml <- "UPSC-0021.Run.xml"
sample.xml <- "UPSC-0021.Sample.xml"
sample.csv <- "ENA_submission_aspwood.csv"
compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0023
#' All samples swapped.
setwd("~/Git/UPSCb/projects/spruce-wood-cross-section/doc/ENA/")
experiment.xml<-"UPSC-0023.Experiment.xml"
run.xml <- "UPSC-0023.Run.xml"
sample.xml <- "UPSC-0023.Sample.xml"
sample.csv <- "ENA_submission_file_info_pabies_full.csv"
U23 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' Fix it
res <- fixXML(experiment.xml,run.xml,sample.csv)
fix.xml<-"UPSC-0023.Experiment-fix.xml"
U23 <- compareXmlToCsv(fix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0025
#' All samples swapped.
setwd("~/Git/UPSCb/projects/spruce-diversity/doc")
i<-25
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "abies.csv"
U25 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' fix - fixed
res <- fixXML(experiment.xml,run.xml,sample.csv)
fix.xml<-paste0("UPSC-00",i,".Experiment-fix.xml")
U25 <- compareXmlToCsv(fix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0026
i<-26
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "wilsonii.csv"
U26<- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0027
#' All samples swapped.
i<-27
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "asperata.csv"
U27 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' fix - fixed
res <- fixXML(experiment.xml,run.xml,sample.csv)
fix.xml<-paste0("UPSC-00",i,".Experiment-fix.xml")
U27 <- compareXmlToCsv(fix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0028
#' All samples swapped.
i<-28
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "crassifolia.csv"
U28 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' fix - fixed
res <- fixXML(experiment.xml,run.xml,sample.csv)
fix.xml<-paste0("UPSC-00",i,".Experiment-fix.xml")
U28 <- compareXmlToCsv(fix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0029
#' All samples swapped.
i<-29
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "engelmannii.csv"
U29 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' fix - fixed
res <- fixXML(experiment.xml,run.xml,sample.csv)
efix.xml<-paste0("UPSC-00",i,".Experiment-fix.xml")
U29 <- compareXmlToCsv(efix.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0030
#' Raw files in different order
i<-30
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "glauca.csv"
U30 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0031
#' Raw files in different order
i<-31
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "glehnii.csv"
U31 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0032
#' Raw files in different order
i<-32
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "jezoensis.csv"
U32 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0033
#' Raw files in different order
i<-33
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "koraiensis.csv"
U33 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0034
#' Raw files in different order
i<-34
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "likiangensis.csv"
U34 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0035
#' Raw files in different order
i<-35
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "mariana.csv"
U35 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0036
#' Raw files in different order
i<-36
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "mexicana.csv"
U36 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0037
#' Raw files in different order
i<-37
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "meyeri.csv"
U37 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0038
#' Raw files in different order
i<-38
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "neoveitchii.csv"
U38 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0039
#' Raw files in different order
i<-39
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "orientalis.csv"
U39 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0040
#' Raw files in different order
i<-40
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "pungens.csv"
U40 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0041
#' Raw files in different order
i<-41
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "retroflexa.csv"
U41 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0042
#' Raw files in different order
i<-42
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "rubens.csv"
U42 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0043
#' Raw files in different order
i<-43
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "schrenkiana.csv"
U43 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0044
#' Raw files in different order
i<-44
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "smithiana.csv"
U44 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0045
#' Raw files in different order
i<-45
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "alcoquiana.csv"
U45 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0046
i<-46
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "aurantiaca.csv"
U46 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0047
#' Raw files in different order
i<-47
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "brachytyla.csv"
U47 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0048
#' Raw files in different order
i<-48
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "breweriana.csv"
U48 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0049
#' Raw files in different order
i<-49
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "chihuahuana.csv"
U49 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0050
#' Raw files in different order
i<-50
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "farreri.csv"
U50 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0051
#' Raw files in different order
i<-51
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "jezoensis-2.csv"
U51 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0052
#' Raw files in different order
i<-52
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "koyamae.csv"
U52 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0053
#' Raw files in different order
i<-53
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "maximowiczii.csv"
U53 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0054
#' to REDO
i<-54
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "obovata.csv"
U54 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0055
i<-55
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "purpurea.csv"
U55 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0056
#' Raw files in different order
i<-56
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "spinulosa.csv"
U56 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0057
#' Raw files in different order
i<-57
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "torano.csv"
U57 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0058
#' Raw files in different order
i<-58
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "wilsonii-2.csv"
U58 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0059
#' Raw files in different order
i<-59
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "abies-2.csv"
U59 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0060
i<-60
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "aurantiaca-2.csv"
U60 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0061
#' Raw files in different order
i<-61
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "alcoquiana-2.csv"
U61 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0062
#' Raw files in different order
i<-62
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "asperata-2.csv"
U62 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0063
#' Raw files in different order
i<-63
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "brachytyla-2.csv"
U63 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0064
#' Raw files in different order
i<-64
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "breweriana-2.csv"
U64 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0065
#' Raw files in different order
i<-65
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "chihuahuana-2.csv"
U65 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0066
#' Raw files in different order
i<-66
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "crassifolia-2.csv"
U66 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0067
#' Raw files in different order
i<-67
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "engelmannii-2.csv"
U67 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0068
#' Raw files in different order
i<-68
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "farreri-2.csv"
U68 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0069
#' Raw files in different order
i<-69
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "glauca-2.csv"
U69 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0070
#' Raw files in different order
i<-70
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "glehnii-2.csv"
U70 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0071
#' Raw files in different order
i<-71
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "jezoensis-3.csv"
U71 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0072
#' Raw files in different order
i<-72
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "koraiensis-2.csv"
U72 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0073
i<-73
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "koyamae-2.csv"
U73 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0074
#' Raw files in different order
i<-74
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "likiangensis-2.csv"
U74 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0075
#' Raw files in different order
i<-75
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "mariana-2.csv"
U75 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0076
#' Raw files in different order
i<-76
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "maximowiczii-2.csv"
U76 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0077
#' Raw files in different order
i<-77
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "mexicana-2.csv"
U77 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0078
#' Raw files in different order
i<-78
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "meyeri-2.csv"
U78 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0079
#' Raw files in different order
i<-79
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "neoveitchii-2.csv"
U79 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0080
i<-80
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "chihuahuana-3.csv"
U80 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0081
i<-81
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "orientalis-2.csv"
U81 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0082
#' Raw files in different order
i<-82
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "pungens-2.csv"
U82 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0083
#' Raw files in different order
i<-83
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "retroflexa-2.csv"
U83 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0084
#' Raw files in different order
i<-84
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "rubens-2.csv"
U84 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0085
#' Raw files in different order
i<-85
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "schrenkiana-2.csv"
U85 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0086
#' Raw files in different order
i<-86
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "smithiana-2.csv"
U86 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0087
i<-87
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "spinulosa-2.csv"
U87 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0088
#' Raw files in different order
i<-88
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "torano-2.csv"
U88 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0089
#' Raw files in different order
i<-89
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "wilsonii-3.csv"
U89 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0090
i<-90
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "omorika.csv"
U90 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0091
i<-91
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "sitchensis.csv"
U91 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0092
i<-92
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "brachytyla-3.csv"
U92 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0093
i<-93
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "obovataHiSeq2000.csv"
U93 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0096
i<-96
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "purpurea-2.csv"
U96 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)

#' ## UPSC-0098
setwd("~/Git/UPSCb/projects/spruce-pine-light/doc/ENA")

i<-98
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "../Spruce-light-project-sample-info.csv"

U98 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
correctRunExperimentIDs(experiment.xml,run.xml,sample.csv)
fix.run.xml <- paste0("UPSC-00",i,".Run-fix.xml")
U98 <- compareXmlToCsv(experiment.xml,fix.run.xml,sample.xml,sample.csv)

#' ## UPSC-0099
i<-99
experiment.xml<-paste0("UPSC-00",i,".Experiment.xml")
run.xml <- paste0("UPSC-00",i,".Run.xml")
sample.xml <- paste0("UPSC-00",i,".Sample.xml")
sample.csv <- "../Pine-light-project-sample-info.csv"
U99 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
correctRunExperimentIDs(experiment.xml,run.xml,sample.csv)
fix.run.xml <- paste0("UPSC-00",i,".Run-fix.xml")
U99 <- compareXmlToCsv(experiment.xml,fix.run.xml,sample.xml,sample.csv)

#' ## UPSC-0101
i<-101
setwd("~/Git/UPSCb/projects/bioimprove/doc/ENA")
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "wt_ena_submission.csv"
U101 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U101

#' ## UPSC-0103
i<-103
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "construct9_ena_submission.csv"
U103 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U103

#' ## UPSC-0104
i<-104
setwd("~/Git/USPCb/projects/arabidopsis-rn4-wgs/doc/ENA")
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Template.csv"
U104 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U104

#' ## UPSC-0105
i<-105
setwd("~/Git/UPSCb/projects/bellini_cambium/doc/ENA")
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Sample-Information_CB_IP-OP42.csv"
correctRunExperimentIDs(experiment.xml,run.xml,sample.csv)
fix.run.xml <- paste0("UPSC-0",i,".Run-fix.xml")
U105 <- compareXmlToCsv(experiment.xml,fix.run.xml,sample.xml,sample.csv)
U105

#' ## UPSC-0106
i<-106
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Sample-Information_CB_IP-T89.csv"
correctRunExperimentIDs(experiment.xml,run.xml,sample.csv)
fix.run.xml <- paste0("UPSC-0",i,".Run-fix.xml")
U106 <- compareXmlToCsv(experiment.xml,fix.run.xml,sample.xml,sample.csv)
U106

#' ## UPSC-0102
i<-102
setwd("~/Git/UPSCb/projects/spruce-flakaliden-metagenomics/doc/ENA")
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Sample_Details.csv"
U102 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile=4)
U102

#' ## UPSC-0108
setwd("~/Git/UPSCb/projects/spruce-suspensor-cells/doc/ENA")
i<-108
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "submission.csv"
U108 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U108

#' ## UPSC-0109
setwd("~/Git/UPSCb/projects/fd-2_LY46_5times_3bioreps/doc/ENA")
i<-109
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Harter_submitting_procedure_BPC6_HiSeq_final.csv"
U109 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile=1)
U109

#' ## UPSC-0110
setwd("~/Git/UPSCb/projects/aba/doc/ENA")
i<-110
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ena_submission.csv"
U110 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile=1)
U110

#' ## UPSC-0111
setwd("~/Git/UPSCb/projects/aspseq/doc/ENA")
i<-111
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Potrs-Illumina-HiSeq-100bp.csv"
U111 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U111

#' ## UPSC-0112
setwd("~/Git/UPSCb/projects/aspseq/doc/ENA")
i<-112
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Illumina-HiSeq-100bp.csv"
U112 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U112

#' ## UPSC-0113
setwd("~/Git/UPSCb/projects/aspseq/doc/ENA")
i<-113
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Illumina-HiSeq-75bp.csv"
U113 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U113

#' ## UPSC-0114
setwd("~/Git/UPSCb/projects/aspseq/doc/ENA")
i<-114
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "454.csv"
U114 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U114

#' ## UPSC-0116
setwd("~/Git/UPSCb/projects/arabidopsis-chip-seq-prn2-hub2/doc/ENA")
i<-116
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "samples.csv"
U116 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile=1)
U116

#' ## UPSC-0117
setwd("~/Git/UPSCb/projects/arabidopsis-porcupine-RNA-Seq/doc/ENA")
i<-117
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "submission.csv"
U117 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U117

#' ## UPSC-0118
setwd("~/Git/UPSCb/projects/arabidopsis-FD/doc/ENA")
i<-118
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "FD_manuscript_reads_submission_RNAseq.csv"
U118 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U118

#' ## UPSC-0119
i<-119
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "FD_manuscript_reads_submission_ChIPseq.csv"
U119 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U119

#' ## UPSC-0121
setwd("~/Git/UPSCb/projects/spombe-pfh1-top1/doc/ENA")
i<-121
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Template for data submission.csv"
U121 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U121

### ONGOING

#' ## UPSC-0122
setwd("~/Documents/Git/UPSCb/projects/spruce-needles/doc/ENA")
i<-122
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.csv <- "UPSC-0122.csv"
correctRunExperimentIDs(experiment.xml,run.xml,sample.csv)
fix.run.xml <- paste0("UPSC-0",i,".Run-fix.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
U122 <- compareXmlToCsv(experiment.xml,fix.run.xml,sample.xml,sample.csv)
U122

#' ## UPSC-0123
i<-123
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "UPSC-0123.csv"
U123 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U123

#' ## UPSC-0125
setwd("~/Git/UPSCb/projects/metagenomics/doc/ENA")
i<-125
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "SUBMISSION_AmpSeqFINAL.csv"
U125 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U125

#' ## UPSC-0126
i<-126
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "SUBMISSION_RNASeqFINAL.csv"
U126 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U126

#' ## UPSC-129
setwd("~/Git/UPSCb/projects/Arabidopsis_RNA_cold_acclimation/doc/ENA")
i<-129
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ENA_submission_Arabidopsis_short-term_cold_stress_data_set.c.csv"
U129 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U129

#' ## UPSC-130
setwd("~/Git/UPSCb/projects/spruce-needles-cold-stress/doc/ENA")
i<-130
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ENA_submission_Spruce_cold_stress_needles_data.csv"
U130 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U130

#' ## UPSC-131
setwd("~/Git/UPSCb/projects/spruce-roots-cold-stress/doc/ENA")
i<-131
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ENA_submission_Spruce_cold_stress_roots_data.csv"
U131 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U131

#' ## UPSC-132
setwd("~/Git/UPSCb/projects/spruce-needles-drought-stress/doc/ENA")
i<-132
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ENA_Submission_PA_Drougth_JH_1.csv"
U132 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U132

#' ## UPSC-134
setwd("~/Git/UPSCb/projects/arabidopsis-aphids/doc/ENA")
i<-134
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "submission.csv"
U134 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U134

#' ## UPSC-135
setwd("~/Git/UPSCb/projects/aspseq/doc/ENA")
i<-135
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "UPSC-0135.csv"
U135 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U135

#' ## UPSC-136
i<-136
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "UPSC-0136.csv"
U136 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U136

#' ## UPSC-137
setwd("~/Git/UPSCb/projects/arabidopsis-serrate-ChIP-Seq/doc/ENA")
i<-137
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "TFL1_manuscript_reads_submission_RNAseq.csv"
U137 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U137

#' ## UPSC-138
i<-138
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "TFL1_manuscript_reads_submission_ChIPseq.csv"
U138 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U138

#' ## UPSC-140
setwd("~/Git/UPSCb/projects/ERF/doc/ENA")
i<-140
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ERF139.csv"
U140 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U140

#' ## UPSC-141
i<-141
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "ERF139SRDX.csv"
U141 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U141

#' ## UPSC-0022 update
setwd("~/Git/UPSCb/projects/arabidopsis-serrate-ChIP-Seq/doc/ENA")
i<-22
experiment.xml<-paste0("UPSC-00",i,"-new.Experiment.xml")
run.xml <- paste0("UPSC-00",i,"-new.Run.xml")
sample.xml <- paste0("UPSC-00",i,"-new.Sample.xml")
sample.csv <- "Corinna_submitting_procedure_SE_new.csv"
U22n <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U22n

#' ## UPSC-146
setwd("~/Git/UPSCb/projects/leaf_development/doc/ENA")
i<-146
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "aspleaf_mRNA.csv"
U146 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U146

#' ## UPSC-147
setwd("~/Git/UPSCb/projects/leaf_development-sRNA/doc/ENA")
i<-147
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Template.csv"
U147 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile=1)
U147

#' ## UPSC-149
setwd("~/Box Sync/UPSC/Facility/Doc/Sara/ENA")
i<-149
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "samples.csv"
U149 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U149

#' ## UPSC-149
setwd("~/Git/UPSCb/projects/ethylene-insensitive/doc/ENA")
i<-151
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Novogene.csv"
U151 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U151

#' ## UPSC-152
i<-152
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Sample-TW-20190418.csv"
U152 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U152

#' ## UPSC-153
setwd("~/Git/UPSCb/projects/arabidopsis-mediator-stress/doc/ENA")
i<-153
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
old.run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "Sample_file.csv2"

r.list <- correctRunExperimentIDs(experiment.xml,old.run.xml,sample.csv,n=1)
run.xml <- paste0("UPSC-0",i,".Run-fix.xml")
U153 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv,nfile = 1)
U153

#' ## UPSC-154
setwd("~/Git/UPSCb/projects/bioimprove/doc/ENA")
i<-154
experiment.xml<-paste0("UPSC-0",i,".Experiment.xml")
run.xml <- paste0("UPSC-0",i,".Run.xml")
sample.xml <- paste0("UPSC-0",i,".Sample.xml")
sample.csv <- "construct12_ena_submission.csv"

U154 <- compareXmlToCsv(experiment.xml,run.xml,sample.xml,sample.csv)
U154



