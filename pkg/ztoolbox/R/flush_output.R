#-------------------------------------------------------------------------
# flush output 
#-------------------------------------------------------------------------
flush_output <- function(dynamic.text,fixed.text){
  message(paste0(dynamic.text,fixed.text),"\r",appendLF=FALSE)
  flush.console()
}


