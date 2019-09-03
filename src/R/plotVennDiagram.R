## helper function
myVenn <- function(sets,plot=TRUE,...){

  ## lib
  library(Vennerable)
  
  ## sanity checks
  stopifnot(is.list(sets))

  stopifnot(length(sets) %in% c(2,3))

  stopifnot(!is.null(names(sets)))
  
  ## the original Venn constructor takes for ever.
  ven <- switch(as.character(length(sets)),
                "2"={
                  common <- intersect(sets[[1]],sets[[2]])
                  set1 <- setdiff(sets[[1]],sets[[2]])
                  set2 <- setdiff(sets[[2]],sets[[1]])
                  ven <- new("Venn",IndicatorWeight=
                             matrix(c(rep(0,3),
                                      1,0,length(set1),
                                      0,1,length(set2),
                                      1,1,length(common)),
                                    byrow=TRUE,
                                    dimnames=list(c("00","10","01","11"),
                                      c(names(sets)[1],names(sets)[2],".Weight")),
                                    ncol=3),
                             IntersectionSets=list(
                               "00"=NULL,
                               "10"=set1,
                               "01"=set2,
                               "11"=common))
                },
                "3"={
                  i12 <- intersect(sets[[1]],sets[[2]])
                  s123 <- intersect(i12,sets[[3]])
                  s1 <- setdiff(sets[[1]],union(sets[[2]],sets[[3]]))
                  s2 <- setdiff(sets[[2]],union(sets[[1]],sets[[3]]))
                  s3 <- setdiff(sets[[3]],union(sets[[1]],sets[[2]]))
                  s12 <- setdiff(i12,s123)
                  s23 <- setdiff(intersect(sets[[2]],sets[[3]]),s123)
                  s13 <- setdiff(intersect(sets[[1]],sets[[3]]),s123)

                  ven <- new("Venn",IndicatorWeight=
                             matrix(c(rep(0,4),
                                      1,0,0,length(s1),
                                      0,1,0,length(s2),
                                      1,1,0,length(s12),
                                      0,0,1,length(s3),
                                      1,0,1,length(s13),
                                      0,1,1,length(s23),
                                      1,1,1,length(s123)),
                                    byrow=TRUE,
                                    dimnames=list(c("000","100","010","110","001","101","011","111"),
                                      c(names(sets)[1],names(sets)[2],names(sets)[3],".Weight")),
                                    ncol=4),
                             IntersectionSets=list(
                               "000"=NULL,
                               "100"=s1,
                               "010"=s2,
                               "110"=s12,
                               "001"=s3,
                               "101"=s13,
                               "011"=s23,
                               "111"=s123))
                })
  if(plot){
    plot(ven,...)
  }
  invisible(ven)
}
