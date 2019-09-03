#' Load the data
load("assembly.rda")

# NOTE:
# ANAM: assembly name
# ADESC: description / title
# StREF: study ref
# SaREF: sample ref
# fFile: faste file
# aFile: agp file
# fmd5: fasta md5
# amd5: agp md5

#' Functions
header <- function(filename){
  write('<?xml version="1.0" encoding="UTF-8"?>\n<ANALYSIS_SET>',
        file=filename)
}

body <- function(filename,DROPDIR,
                 ANAM,ADESC,StREF,
                 SaREF,fFile,fmd5,cFile,cmd5){
  write(paste(
    paste0('<ANALYSIS alias="',ANAM,'" center_name="UPSC">'),
    paste0("<TITLE>",ADESC,"</TITLE>"),
    paste0("<DESCRIPTION>",ADESC,"</DESCRIPTION>"),
    paste0('<STUDY_REF refname="',StREF,'" refcenter="UPSC"/>'),
    paste0('<SAMPLE_REF refname="',SaREF,'" refcenter="UPSC"/>'),
    "<ANALYSIS_TYPE>",
    "<SEQUENCE_ASSEMBLY>",
    paste0("<NAME>",ANAM,"</NAME>"),
    "<PARTIAL>0</PARTIAL>",
    "<COVERAGE>50</COVERAGE>",
    "<PROGRAM>SPades v. 3.7.0</PROGRAM>",
    "<PLATFORM>ILLUMINA</PLATFORM>",
    "</SEQUENCE_ASSEMBLY>",
    "</ANALYSIS_TYPE>",
    "<FILES>",
    paste0('<FILE filename="',file.path(DROPDIR,fFile),
           '" filetype="chromosome_fasta" checksum_method="MD5" checksum="',
           fmd5, '" />'),
    paste0('<FILE filename="',file.path(DROPDIR,cFile),
           '" filetype="chromosome_list" checksum_method="MD5" checksum="',
           cmd5,'" />'),
    "</FILES>",
    "</ANALYSIS>",
    sep ="\n"),
    file=filename,
    append = TRUE)
}

footer <- function(filename){
  write('</ANALYSIS_SET>',
        file=filename,
        append = TRUE)
}

#' Process
filename="UPSC-0094.Analysis.xml"
header(filename)
dev.null <- apply(df,1,function(d){
  body(filename,"UPSC-0094",
       sub("\\.fasta\\.gz","",d["fFile"]),
       paste("Chloroplast assembly from a",d["SNAM"]),
       d["StRef"],d["SaREF"],d["fFile"],d["fMd5"],
       d["cFile"],d["cMd5"]
       )
})
footer(filename)
