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

##' 
##' Reads in the JHUCSSE total case count data up
##' until (and including) a given dat.
##' 
##' @param date the last day to consider data from
##'
##' @return a data frame with the basic data.
##' 
read_JHUCSSE_cases <- function(date) {

  ## first get a list of all of the files in the directory
  ## starting with "JHUCSSE Total Cases"
  file_list <- list.files("data","JHUCSSE Total Cases",
                          full.names = TRUE)
  
  file_list <- rev(file_list)
  print(file_list)
  
  ##Now combine them into one data frame
  rc <- NULL
  
  for (file in file_list) {
    tmp <- read_csv(file)%>%
      rename(Province_State=`Province/State`)%>%
      rename(Update = `Last Update`) %>%
      mutate(Update=lubridate::parse_date_time(Update, c("%m/%d/%Y %I:%M %p", "%m/%d/%Y %H:%M", "%m/%d/%y %I:%M %p")))
    
    if("Country"%in%colnames(tmp)) {
      tmp <- rename(tmp, Country_Region=Country)
    } else {
      tmp <- rename(tmp, Country_Region=`Country/Region`)
    }
     
    rc <-bind_rows(rc,tmp)
  }
  
  ##Now drop any after the date given
  rc <- rc%>%filter(Update<=date)
  
  return(rc)
}