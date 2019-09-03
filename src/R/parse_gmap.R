#'---
#'title: Simple parsing utility for standard gmap summary output
#'author: Bastian Schiffthaler
#'---
#'
#'This function utilizes a Python based lexic parser, which matches regular expression in the 
#'standard output of gmap and parses them into a GRanges object with most of the info as 
#'additional columns. 
parse.gmap <- function(gmapfile){
  require(rPython)
  require(GenomicRanges)
  tic <- Sys.time()
  #The Python code is currently in a Python source file but could be hardcoded here as a 
  #character object to make this script self-contained
  parser.py <- "/mnt/picea/home/bastian/Git/UPSCb/src/python/parse_gmap.py"
  #The gmap input file is assignedas a variable in the Python buffer
  python.assign('gmapfile',gmapfile)
  #The file is read as character and executed in the Python environment
  message("Strting to read input file:",gmapfile)
  python.exec(readChar(parser.py,file.info(parser.py)$size))
  message("Done reading file... converting to GRanges")
  #We can retrieve the resulting object from the Python buffer and get a nested list in return
  gmap_dict <- python.get(var.name = "DICT")
  #Getting transcript matches
  paths <- unlist(lapply(gmap_dict,length))
  transcripts <- rep(names(paths),paths)
  #We unlist one level and rbind the list to a char matrix into as.data.frame
  gmap_dict <- unlist(gmap_dict,recursive = FALSE)
  gmap_dict <- as.data.frame(do.call(rbind, gmap_dict), stringsAsFactors=FALSE)
  #Since IRanges will not allow ranges from a larger to a smaller number, we can logically return
  #the strand and later swap in the GRanges creator based on which number is bigger/smaller
  strand <- ifelse(as.integer(gmap_dict$Genome.start)>as.integer(gmap_dict$Genome.end),"-","+")
  gr <- GRanges(seqnames = gmap_dict$Genome.loc,ranges = IRanges(
    start = ifelse(as.integer(gmap_dict$Genome.start)>as.integer(gmap_dict$Genome.end),
                   as.integer(gmap_dict$Genome.end),as.integer(gmap_dict$Genome.start)),
    end = ifelse(as.integer(gmap_dict$Genome.end)<as.integer(gmap_dict$Genome.start),
                 as.integer(gmap_dict$Genome.start),as.integer(gmap_dict$Genome.end))
    ), strand = strand, "Transcript.match" = transcripts, 
    "Transcript.match.start" = as.integer(gmap_dict$Transcript.start),
    "Transcript.match.end" = as.integer(gmap_dict$Transcript.start),
    "Coverage" = as.numeric(gmap_dict$Coverage),
    "Trimmed.coverage" = as.numeric(gmap_dict$Trimmed.coverage),
    "Percent.identity" = as.numeric(gmap_dict$Percent.identity),
    "Matches" = as.integer(gmap_dict$Matches),
    "Mismatches" = as.integer(gmap_dict$Mismatches),
    "Indels" = as.integer(gmap_dict$Indels)
    )
  message("Done in ", round(difftime(Sys.time(), tic, units = "secs"),2)," seconds.")
  return(gr)
}

