##' Function to read in and do some minor cleaning on the
##' "Kudos" data.
##' 
##' @param filename the file name
##' 
##' @return a data frame with the basic data.
##' 
readKudos2 <- function (filename) {
  require(tidyverse)
  rc <- read_csv(filename, col_types = 
                   cols(date=col_date("%m/%d/%Y"),
                        gender=col_factor(),
                        symptom_onset = col_date("%m/%d/%Y")))
  
 
  return(rc)
}

##' Same as above but for older Kudos file layout.
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