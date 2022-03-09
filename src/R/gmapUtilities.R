#' # Setup
#' Libraries
suppressPackageStartupMessages({
  require(dplyr)
  require(GenomicRanges)
  require(igraph)
  require(magrittr)
  require(readr)
  require(tidyr)
})

#' # Class
setClass("GmapGff3",
         contains="tbl_df")

#' # Generics
#' ## S4 Export-able
setGeneric(name="read_gmap_gff3",
           def=function(file=character(0),
                        type=c("gene","mRNA","exon"),
                        guess_type=logical(0)){
             standardGeneric("read_gmap_gff3")})

#' ## S4 Internal (workers)
setGeneric(name=".read_gmap_gff3_gene",
           def=function(file=character(0),...){
             standardGeneric(".read_gmap_gff3_gene")})

setGeneric(name=".read_gmap_gff3_mrna",
           def=function(file=character(0),...){
             standardGeneric(".read_gmap_gff3_mrna")})

setGeneric(name=".read_gmap_gff3_exon",
           def=function(file=character(0),...){
             standardGeneric(".read_gmap_gff3_exon")})

setGeneric(name=".basic_stats",
           def=function(obj=NULL){
             standardGeneric(".basic_stats")})

#' # Methods
#' ## S4 Export-able
setMethod(f = "read_gmap_gff3", 
          signature = c("character"),
          definition = function(
            file=character(0),
            type=c("gene","mRNA","exon"),
            guess_type=TRUE
          ){
            
            # first sanity
            stopifnot(file.exists(file))
            
            ## dispatch, also sanitize
            as(switch(match.arg(type),
                   "gene" = .read_gmap_gff3_gene(file,guess_type),
                   "mRNA" = .read_gmap_gff3_mrna(file,guess_type),
                   "exon" = .read_gmap_gff3_exon(file,guess_type)),
               "GmapGff3")
          })

#' ## S4 coercion
setAs("GmapGff3","GRanges",
      function(from){
        return(GRanges(seqnames=from$Chr,
                ranges=IRanges(start=from$Start,end=from$End),
                strand=from$Strand,ID=from$ID))
      })

#' ## S4 print
setMethod(f="print",
          signature="GmapGff3",
          definition=function(x){
            print(as(x,"tbl_df"))
            cat("# Basic stats\n")
            cat(.basic_stats(x),"\n")
          })

#' ## S4 width
setMethod(f="width",
          signature="GmapGff3",
          definition=function(x){
            len <- x$End-x$Start+1
            if(!is.null(x$Type)){
              if(any(x$Type != "uniq")){
                len <- sapply(split(len,x$ID),sum)
              }
            }
            return(len)
          })

#' ## S4 internal
setMethod(f = ".read_gmap_gff3_gene", 
          signature = c("character"),
          definition = function(
            file=character(0),
            ...
          ){
            return(.read(file) %>%
                     filter(Type=="gene") %>% 
                     .select_gene() %>% 
                     .guess_type(...))
          })

setMethod(f = ".read_gmap_gff3_mrna", 
          signature = c("character"),
          definition = function(
            file=character(0),
            ...
          ){
            return(.read(file) %>%
                     filter(Type=="mRNA") %>% 
                     .select_mrna() %>% 
                     .guess_type(...))
          })

setMethod(f=".basic_stats",
          signature="ANY",
          definition=function(obj){
            stats <- .bp_stats(width(obj))
            sprintf("They cover %s +/- %s bp on average at the genomic level (min=%s, max=%s)",
                    stats[1],stats[2],stats[3],stats[4])
          })

#' ## Internal functions
".read" <- function(file){
  # file will already be sanitized for existence
  return(read_tsv(file,
                  col_names=c("Chr","Source","Type","Start",
                              "End","Score","Strand","Frame","Attrs"),
                  comment="#",
                  show_col_types = FALSE))
}

".select_gene" <- function(tbl){
  return(tbl %>% mutate(Chr=as.factor(Chr),
                        Strand=as.factor(Strand),
                        ID=gsub("ID=|\\.path\\d+;.*","",Attrs)) %>% 
           select(Chr,Start,End,Strand,ID))
}

".select_mrna" <- function(tbl){
  return(tbl %>% mutate(Chr=as.factor(Chr),
                       Strand=as.factor(Strand))
         %>% separate(Attrs,into=c(NA,"ID",NA,"Name",NA,"Parent",
                                   NA,"Dir",NA,"Coverage",NA,"Identity",
                                   NA,"Matches",NA,"Mismatches",NA,"Indels",NA,NA),
                      sep="=|;") %>%  
           mutate(Coverage=parse_double(Coverage),
                  Identity=parse_double(Identity),
                  Matches=parse_integer(Matches),
                  Mismatches=parse_integer(Mismatches),
                  Indels=parse_integer(Indels)) %>% 
           select(-Source,-Type,-Score,-Frame))
}

".guess_type" <- function(tbl,guess_type){
  if(guess_type){
    return(left_join(tbl,
      tbl %>% select(ID) %>% 
        group_by(ID) %>% 
        dplyr::count() %>% 
        mutate(Type=ifelse(n==1,"uniq",ifelse(n>2,"mult","multOrTransloc"))) %>% 
        select(ID,Type),by="ID"))
  } else {
    return(tbl)
  }
}

".bp_stats" <- function(len){
  return(c(round(mean(len)),
                round(sd(len)),
                min(len),
                max(len)))
}
  
  
  