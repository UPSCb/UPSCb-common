##########
### Dynamically extend the PATH in R to access modules via system()
##########

module.load <- function(args=""){
  require(stringr)
  
  # If no MODULEPATH is set, assume default
  if(Sys.getenv("MODULEPATH")==""){
  modulefile <- readLines("/mnt/picea/storage/Modules/default/init/.modulespath")
  modulepath <- sapply(strsplit(str_extract(modulefile,"^\\/.*\\t")[!is.na(str_extract(modulefile,"^\\/.*\\t"))],split="\t"),"[[",1)
  modulepath <- paste(modulepath,collapse=":")
  Sys.setenv("MODULEPATH" = modulepath)
  }
  # Capture response of the module command
  system.return <- system(paste("/mnt/picea/Modules/default/bin/modulecmd sh load",args), intern = TRUE)
  
  # If the MODULEPATH is to be extended, detect the system return for the string:
  if(str_detect(system.return,"export MODULEPATH")){
  system.return <- str_replace(str_extract(system.return,";MODULEPATH=.*?;"),";MODULEPATH=","")
  system.return <- str_replace(system.return,";","")
  Sys.setenv("MODULEPATH" = paste(Sys.getenv("MODULEPATH"),system.return,sep=":"))
  # Make sure we remove dupes
  long.path <- Sys.getenv("MODULEPATH")
  short.path <- paste(unique(str_replace(unlist(strsplit(long.path,":"))," +","")),collapse=":")
  Sys.unsetenv("MODULEPATH")
  Sys.setenv("MODULEPATH" = short.path)
  
  } else if (str_detect(system.return,"export PATH")) {
    # Otherwise PATH should be extended
    system.return <- str_replace(str_extract(system.return,";PATH=.*?;"),";PATH=","")
    system.return <- str_replace(system.return,";","")
    Sys.setenv("PATH" = paste(Sys.getenv("PATH"),system.return,sep=":"))
    # Make sure we remove dupes
    long.path <- Sys.getenv("PATH")
    short.path <- paste(unique(str_replace(unlist(strsplit(long.path,":"))," +","")),collapse=":")
    Sys.unsetenv("PATH")
    Sys.setenv("PATH" = short.path)
  
    } else {
    stop("Don't know how to process that: no PATH or MODULEPATH directive detected")
  }
};
