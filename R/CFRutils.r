#' A function to gather data into right format for model
#' ##' @param date todays date
##'
##' @return a list with the data to fit CFR model
##' @example  mod2_nohubei <- moddata(ISOdate(2020,2,9, hour = 23, min = 59, tz = "EST") , withHubei = FALSE)
##' 
##' 
##' 
`%nin%` = Negate(`%in%`)

moddata <- function(date, withHubei = TRUE ){ # outsideChina = FALSE, ChinaOnly = FALSE, , HubeiOnly = FALSE
  # 
  # if(outsideChina == TRUE){
  #   ChinaOnly <- FALSE
  #   withHubei <- FALSE
  #   HubeiOnly <- FALSE
  # }
  #read in JHUCSSE data
  jhucsse <- read_JHUCSSE_cases(date, 
                                append_wiki = TRUE) 
  
  ##Continue to filter to China for the moment (we can turn this off)

  jhucsse_use <- jhucsse %>% 
    filter(Country_Region%in%c("Mainland China", "Macau", "Hong Kong")) 
  
 #jhucsse_use <- jhucsse
 # 
 # if(outsideChina == TRUE){
 #   jhucsse_use <- jhucsse_use %>% 
 #     filter(Country_Region%nin%c("Mainland China", "Macau", "Hong Kong")) 
 # }
 # 
 # if(ChinaOnly == TRUE){
 #   jhucsse_use <- jhucsse_use %>% 
 #      filter(Country_Region%in%c("Mainland China", "Macau", "Hong Kong")) 
 # }
 #  
 #  if(withHubei == FALSE){
 #    jhucsse_use <- jhucsse_use %>% filter(Province_State != "Hubei")
 #  }
 #  
 #  if(HubeiOnly == TRUE){
 #    jhucsse_use <- jhucsse %>% filter(Province_State == "Hubei")
 #  }
  
  
 jhucsse_use$Province_State <- as.factor(jhucsse_use$Province_State)
 jhucsse_use$loctest <- as.numeric(jhucsse_use$Province_State)
  #estimate daily incidence
  incidence_data <- est_daily_incidence(jhucsse_use,
                                        ISOdate(2019,12,1),
                                        date)
  
  incidence_data<-
    incidence_data%>%filter(Date>"2020-1-1")
 
 
  #Add columns with indices for time and location
  incidence_data$loc <-
    as.numeric(incidence_data$Province_State)
  incidence_data$t <- as.numeric(as.Date(incidence_data$Date))-
    min(as.numeric(as.Date(incidence_data$Date))) + 1 
  
  #Max day of follow-up
  Tmax <- max(incidence_data$t)
  #number of locations considered
 # L <- length(unique(incidence_data$loc))#
  L <- max(incidence_data$loc)
  
  cases <- matrix(0,nrow=Tmax, ncol=L)
  
  for (i in 1:nrow(incidence_data)) {
    cases[incidence_data$t[i], incidence_data$loc[i]] <-
      incidence_data$Incidence[i]
  }
  
  #estimate daily recoveries
  recov_data <- est_daily_recovered(jhucsse_use,
                                    ISOdate(2019,12,1),
                                    date,
                                    na_to_zeros = TRUE)
  
  recov_data<-
    recov_data%>%filter(Date>"2020-1-1") %>%
    drop_na(Recovered)

  #Add columns with indices for time and location
  recov_data$loc <-
    as.numeric(recov_data$Province_State)
  recov_data$t <- as.numeric(as.Date(recov_data$Date))-
    min(as.numeric(as.Date(recov_data$Date))) + 1
  
  
  recovered <- matrix(0,nrow=Tmax, ncol=L)
  for (i in 1:nrow(recov_data)) {
    recovered[recov_data$t[i], recov_data$loc[i]] <-
       recov_data$Recovered[i]
  }
  
  
  #estimate daily deaths
  death_data <- est_daily_deaths(jhucsse_use,
                                 ISOdate(2019,12,1),
                                 date)
  
  death_data<-
    death_data%>%filter(Date>"2020-1-1") %>% 
    drop_na(Deaths)
  
  
  #Add columns with indices for time and location
  death_data$loc <-
    as.numeric(as.numeric(death_data$Province_State))
  death_data$t <- as.numeric(as.Date(death_data$Date))-
    min(as.numeric(as.Date(incidence_data$Date))) + 1
  
  
  deaths <- matrix(0,nrow=Tmax, ncol=L)
  
  for (i in 1:nrow(death_data)) {
    deaths[death_data$t[i], death_data$loc[i]] <-
      death_data$Deaths[i]
  }
  
  #housekeeping: ensure nothing is negative
  deaths[deaths<0] <- 0
  cases[cases <= 0] <- 0.0001
  recovered[recovered<0] <- 0
  
  #make indicator of whether or not position x is hubei
  w <- vector(length = ncol(cases))
  w[unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))] <- 1
  # 
  #omit hubei if desired
  if(withHubei == FALSE){
    deaths <- deaths[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
    recovered <- recovered[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
    cases <- cases[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
    w <- w[ -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
    L <- length(w) #remake L if omitting hubei
  }

  V <- nrow(cases) #max infection duration
  
  #gather data into a list
  cfrmdl_dat <- list(T=Tmax, L=L,
                     c = cases,
                     d = round(deaths),
                     r = round(recovered),
                     w = array(w), 
                     V = V)
  
  return(cfrmdl_dat)
}