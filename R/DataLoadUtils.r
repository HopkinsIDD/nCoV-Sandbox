##' Function to read in and do some minor cleaning on the
##' "Kudos" data.
##' 
##' @param filename the file name
##' 
##' @return a data frame with the basic data.
##' 
readKudos <- function (filename) {
  require(tidyverse)
  rc <- read_csv(filename, col_types = 
                   cols(Date=col_date("%m/%d/%Y"),
                        Gender=col_factor(),
                        `Symptom onset (approximate)` = col_date("%m/%d/%Y")))
  
  #extract death data
  rc <- mutate(rc, dead=str_detect(Summary,"death")) %>%
    rename(onset=`Symptom onset (approximate)`)
  
  return(rc)
}