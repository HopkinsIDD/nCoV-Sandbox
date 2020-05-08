#' A function to gather data into right format for model
#' ##' @param date todays date
##'
##' @return a list with the data to fit CFR model
##' @example  mod2_nohubei <- moddata(ISOdate(2020,2,9, hour = 23, min = 59, tz = "EST") , withHubei = FALSE)
##' 
##' 
##' 
`%nin%` = Negate(`%in%`)

moddata <- function(date, withHubei = TRUE){ # outsideChina = FALSE, ChinaOnly = FALSE, , HubeiOnly = FALSE
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
  jhucsse_use <- jhucsse
  
   # jhucsse_use <- jhucsse %>%
   #                  group_by(Country_Region) %>% 
   #                  summarise(tConfirmed = sum(Confirmed), tDeaths = sum(Deaths), tRecovered = sum(Recovered))
   # 
   # jhucsse_use <- jhucsse_use %>% 
   #   rename(Confirmed=tConfirmed, Deaths=tDeaths, Recovered=tRecovered)
  
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


moddata_country <- function(date){ # outsideChina = FALSE, ChinaOnly = FALSE, , HubeiOnly = FALSE
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
 
  
  jhucsse_use <- jhucsse 
  #clean up some countries
  jhucsse_use$Country_Region <- ifelse(jhucsse_use$Country_Region=="France"&jhucsse_use$Province_State!="France", jhucsse_use$Province_State, jhucsse_use$Country_Region)
  #separrrate out hubei
  jhucsse_use$Country_Region <- ifelse(jhucsse_use$Province_State=="Hubei", jhucsse_use$Province_State, jhucsse_use$Country_Region)
 
  jhucsse_use <- jhucsse_use %>% 
    mutate(Country_Region=replace(Country_Region, Country_Region=="United States", "US")) %>%
    mutate(Country_Region=replace(Country_Region, Country_Region=="United Kingdom", "UK")) %>%
    mutate(Country_Region=replace(Country_Region, Country_Region=="China", "Mainland China")) 
  
 # jhucsse_use$Country_Region <- ifelse(jhucsse_use$Province_State=="Hubei", jhucsse_use$Province_State, jhucsse_use$Country_Region)
  
  # %>%
  #                  group_by(Country_Region, Update) %>%
  #                  summarise(tConfirmed = sum(Confirmed), tDeaths = sum(Deaths), tRecovered = sum(Recovered))

  # jhucsse_use <- jhucsse_use %>%
  #   rename(Confirmed=tConfirmed, Deaths=tDeaths, Recovered=tRecovered, Province_State = Country_Region) #abusing the naming convention here to be able to use est_daily_incidence function
  # 
  # jhucsse_use <- jhucsse_use[!is.na(jhucsse_use$Province_State),]
  # jhucsse_use$Province_State <- ifelse(jhucsse_use$Province_State=="Viet Nam", "Vietnam", jhucsse_use$Province_State)
  
  jhucsse_use$Province_State <- as.factor(jhucsse_use$Province_State)
  jhucsse_use$loctest <- as.numeric(jhucsse_use$Province_State)
  #estimate daily incidence
 
  
  
  # if(smooth ==FALSE){
  #   tmp_dt_seq <- seq( ISOdate(2020,1,1), date, "days")
  #   analyze <-   jhucsse_use %>% replace(is.na(.), 0)
  #   incidence_data <- analyze %>% 
  #     arrange(Country_Region, Province_State, Admin2, Update) %>% 
  #     group_by(Country_Region, Province_State, Admin2, Update) %>%
  #     summarise(cumcases = sum(Confirmed)) %>%
  # #     ungroup(Country_Region, Update) %>% 
  # #     arrange(Country_Region, Update) %>% 
  # #     group_by(Country_Region) %>% 
  #      mutate(cs = cumcases - lag(cumcases, n=1, default = 0)) %>% 
  #  ungroup(Country_Region, Province_State, Admin2, Update) %>%
  #    group_by(Country_Region, Update) %>% 
  #    summarize(c = sum(cs))
  #     
  #     
  # }
  

    incidence_data <- est_daily_incidence_country(jhucsse_use,
                                                  ISOdate(2019,12,1),
                                                  date, 
                                                  na_to_zeros = FALSE)
    

  
  # else{
  #   incidence_data <- 
  # }
  
  incidence_data<-
    incidence_data%>%filter(Date>"2020-1-1")
  
  incidence_data_country <- incidence_data %>% 
    group_by(Country_Region, Date) %>% 
    summarise(tIncidence = sum(Incidence)) %>% 
    rename(Incidence=tIncidence)
  
  
  #Add columns with indices for time and location
  incidence_data_country$loc <-
    as.numeric(factor(incidence_data_country$Country_Region))
  incidence_data_country$t <- as.numeric(as.Date(incidence_data_country$Date))-
    min(as.numeric(as.Date(incidence_data_country$Date))) + 1 
  
  #Max day of follow-up
  Tmax <- max(incidence_data_country$t)
  #number of locations considered
  # L <- length(unique(incidence_data$loc))#
 L <- max(incidence_data_country$loc, na.rm = TRUE)
  
  cases <- matrix(0,nrow=Tmax, ncol=L)
  
  for (i in 1:nrow(incidence_data_country)) {
    cases[incidence_data_country$t[i], incidence_data_country$loc[i]] <-
      incidence_data_country$Incidence[i]
  }
  
  #estimate daily recoveries
  recov_data <- est_daily_recovered_country(jhucsse_use,
                                    ISOdate(2019,12,1),
                                    date,
                                    na_to_zeros = FALSE)
  
  recov_data<-
    recov_data%>%filter(Date>"2020-1-1") %>%
    drop_na(Recovered)
  
  recov_data_country <- recov_data %>% 
    group_by(Country_Region, Date) %>% 
    summarise(tRecovered = sum(Recovered)) %>% 
    rename(Recovered=tRecovered)
  
  
  
  #Add columns with indices for time and location
  recov_data_country$loc <-
    as.numeric(factor(recov_data_country$Country_Region))
  recov_data_country$t <- as.numeric(as.Date(recov_data_country$Date))-
    min(as.numeric(as.Date(incidence_data_country$Date))) + 1
  
  
  recovered <- matrix(0,nrow=Tmax, ncol=L)
  for (i in 1:nrow(recov_data_country)) {
    recovered[recov_data_country$t[i], recov_data_country$loc[i]] <-
      recov_data_country$Recovered[i]
  }
  
  
  #estimate daily deaths
  death_data <- est_daily_deaths_country(jhucsse_use,
                                 ISOdate(2019,12,1),
                                 date)
  
  death_data<-
    death_data%>%filter(Date>"2020-1-1") %>% 
    drop_na(Deaths)
  
  death_data_country <- death_data %>% 
    group_by(Country_Region, Date) %>% 
    summarise(tDeaths = sum(Deaths)) %>% 
    rename(Deaths=tDeaths)
  
  #Add columns with indices for time and location
  death_data_country$loc <-
    as.numeric(as.factor(death_data_country$Country_Region))
  death_data_country$t <- as.numeric(as.Date(death_data_country$Date))-
    min(as.numeric(as.Date(incidence_data_country$Date))) + 1
  
  
  deaths <- matrix(0,nrow=Tmax, ncol=L)
  
  for (i in 1:nrow(death_data_country)) {
    deaths[death_data_country$t[i], death_data_country$loc[i]] <-
      death_data_country$Deaths[i]
  }
  
  #housekeeping: ensure nothing is negative
  deaths[deaths<0] <- 0
  cases[cases <= 0] <- 0.0001
  recovered[recovered<0] <- 0
  
  #make indicator of whether or not position x is hubei
 # w <- vector(length = ncol(cases))
  w <- rep(1, ncol(cases))
  w[unlist(c(head(incidence_data_country[incidence_data_country$Country_Region == "South Korea",  "loc"], n=1), 
             head(incidence_data_country[incidence_data_country$Country_Region == "Diamond Princess",  "loc"], n=1),
             head(incidence_data_country[incidence_data_country$Country_Region == "Mainland China",  "loc"], n=1)))] <- 0
  # 
  #omit hubei if desired
  # if(withHubei == FALSE){
  #   deaths <- deaths[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
  #   recovered <- recovered[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
  #   cases <- cases[, -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
  #   w <- w[ -unlist(head(incidence_data[incidence_data$Province_State == "Hubei",  "loc"], n=1))]
  #   L <- length(w) #remake L if omitting hubei
  # }
  # 
  V <- nrow(cases) #max infection duration
  
  #gather data into a list
  cfrmdl_dat <- list(T=Tmax, L=L,
                     c = cases,
                     d = round(deaths),
                     r = round(recovered),
                     w = array(w), 
                     V = V)
  
  #make a lookup table
  locs <- unique(incidence_data_country$loc)
  names(locs) <- unique(incidence_data_country$Country_Region)
  locs
  return(list(cfrmdl_dat, locs))
}





### Likelihood functions for MLE approach ###

##' likelihood function for exponential dist with only aggregate data
##'
##' @param params startinv values for parameters
##' @param data dataset including confirmed cases, deaths, recovereds on each day, T (# time points), L (# areas), td (vector of individual death times), and tr (vector of individual recovery times)
##'        in epi.curve with missing observations.
##' @return -1 * likelihood
LLexp.agg <- function(params, data){
  # initialization parameters
  h1 <- exp(params[1])
  h2 <- exp(params[2])
  
  # format input data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  
  # initialize expected values
  expected_riskset = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0's
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_deaths to all 0's
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_recovereds to all 0's
  
  # compute expected values
  for (j in 1:L) {  # Loops through countries (L)
    expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
    expected_deaths[1,j] = expected_riskset[1,j] * h1 
    expected_recovereds[1,j] = expected_riskset[1,j] * h2 
    for (t in 2:T) {  # Loops through times within each country
      expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                 (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
      expected_deaths[t,j] = expected_riskset[t,j] *  h1  
      expected_recovereds[t,j] = expected_riskset[t,j] *  h2 
    } 
  }
  
  # likelihood formula
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + 
              (-expected_recovereds + r*log(expected_recovereds)))
  
  # returns negative since optim() minimizes (rather than maximizes)
  return(-LL)
}


##' likelihood function for exponential dist with aggregate + individual level data
##'
##' @param params startinv values for parameters
##' @param data dataset including confirmed cases, deaths, recovereds on each day, T (# time points), L (# areas), td (vector of individual death times), and tr (vector of individual recovery times)
##'        in epi.curve with missing observations.
##'
##' @return likelihood
LLexp.ind <- function(params, data){
  # initialization parameters
  h1 <- exp(params[1])
  h2 <- exp(params[2])
  
  # format input data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  td <- data$td  # individual-level data
  tr <- data$tr  # individual-level data
  
  # initialize expected values
  expected_riskset = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  
  #compute expected values for aggregate data
  for (j in 1:L) {  # Loops through countries
    expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
    expected_deaths[1,j] = expected_riskset[1,j] * h1 
    expected_recovereds[1,j] = expected_riskset[1,j] * h2 
    for (t in 2:T) {  # Loops through all time points
      expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                 (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
      expected_deaths[t,j] = expected_riskset[t,j] *  h1  
      expected_recovereds[t,j] = expected_riskset[t,j] *  h2 
    } 
  }
  
  # compute likelihood contributions for individual data
  lftd <- log(h1*exp(-(h1+h2)*td))
  lftr <- log(h2*exp(-(h1+h2)*tr))
  
  # likelihood formula
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(lftd) + sum(lftr)
  
  # returns negative since optim() minimizes (rather than maximizes)
  return(-LL)
}



##' likelihood function for exponential dist with aggregate + individual level data with random effect (NEEDS WORK)
##'
##' @param params startinv values for parameters
##' @param data dataset including confirmed cases, deaths, recovereds on each day, T (# time points), L (# areas), 
##' td (vector of individual death times), and tr (vector of individual recovery times), wdi (vector of locations of indiidual deaths), wri (vector of locations of individual recoveries)
##'
##' @return likelihood

LLexp.ind.random <- function(params, data){
  #parameters
  loglambda1i <- c(params[3:(2+L)])
  loglambda2i <- c(params[(3+L):(2+2*L)])
  loglambda1 <- params[1]
  loglambda2 <- params[2]
  
  # data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  td <- dat$td
  tr <- dat$tr
  wdi <- dat$wdi
  wdr <- dat$wdr
  rp <- 0 # how much to penalize towards average values of loglambda1 and 2 (commented out to let rp be estimated)
  
  # initialize expected values
  expected_riskset = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  
  #compute h1 and h2
  h1 = exp(loglambda1 + loglambda1i)
  h2 = exp(loglambda2 + loglambda2i)
  
  #compute expected values for aggregate data
  for (j in 1:L) {
    
    expected_riskset[1,j] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
    expected_deaths[1,j] = expected_riskset[1,j] * h1[j] 
    expected_recovereds[1,j] = expected_riskset[1,j] * h2[j] 
    
    for (t in 2:T) {
      expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - (expected_riskset[t-1,j]*h1[j] +  (expected_riskset[t-1,j]* h2[j])))
      expected_deaths[t,j] = expected_riskset[t,j] *  h1[j]  
      expected_recovereds[t,j] = expected_riskset[t,j] *  h2[j] 
    } 
  }
  
  #compute likelihood contributions for individual data
  lftd <- log(h1[wdi]*exp(-(h1[wdi]+h2[wdi])*td))
  lftr <- log(h2[wri]*exp(-(h1[wri]+h2[wri])*tr))
  
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(lftd) + sum(lftr) 
  - sum(.5*rp*(loglambda1i - 0)^2) - sum(.5*rp*(loglambda2i-0)^2)
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}

##' likelihood function for weibull dist with aggregate + individual level data
##'
##' @param params starting values for parameters
##' @param data dataset including confirmed cases, deaths, recovereds on each day, T (# time points), L (# areas), td (vector of individual death times), and tr (vector of individual recovery times)
##'        in epi.curve with missing observations.
##'
##' @return likelihood


LLwbl.agg <- function(params, data){
  #parameters
  loglambda1 <- params[1]
  loglambda2 <- params[2]
  logalpha1 <- params[3]
  logalpha2 <- params[4]
  
  # data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  td <- dat$td
  tr <- dat$tr
  
  # initialize expected values
  
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  alpha1 <- exp(logalpha1)
  alpha2 <- exp(logalpha2)
  
  k <- c(1:T)
  #compute h1 and h2, matarices of hazard functions by day since infection
  h1 = alpha1*exp(loglambda1)*k^(alpha1-1)
  h2 = alpha2*exp(loglambda2)*k^(alpha2-1)
  
  #compute expected values for aggregate data
  for (j in 1:L) {
    expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
    expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
    expected_deaths[1,j] = expected_riskset[1,1] * h1[1] 
    expected_recovereds[1,j] = expected_riskset[1,1] * h2[1]
    
    for (t in 2:T) {
      expected_riskset[1,t] = c[t,j] 
      expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] 
      expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1]
      
      for(k in 2:t){
        
        expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
        expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
        expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k]);
        
      }
    }
  }
  
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds)))
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}


LLwbl.ind <- function(params, data){
  #parameters
  loglambda1 <- params[1]
  loglambda2 <- params[2]
  logalpha1 <- params[3]
  logalpha2 <- params[4]
  
  # data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  td <- dat$td
  tr <- dat$tr
  
  # initialize expected values
  
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  alpha1 <- exp(logalpha1)
  alpha2 <- exp(logalpha2)
  
  k <- c(1:T)
  #compute h1 and h2, matarices of hazard functions by day since infection
  h1 = alpha1*exp(loglambda1)*k^(alpha1-1)
  h2 = alpha2*exp(loglambda2)*k^(alpha2-1)
  
  #compute expected values for aggregate data
  for (j in 1:L) {
    expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
    expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
    expected_deaths[1,j] = expected_riskset[1,1] * h1[1] 
    expected_recovereds[1,j] = expected_riskset[1,1] * h2[1]
    
    for (t in 2:T) {
      expected_riskset[1,t] = c[t,j] 
      expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] 
      expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1]
      
      for(k in 2:t){
        
        expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
        expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
        expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k]);
        
      }
    }
  }
  
  #compute likelihood contributionss for individual data
  lftd <- log(alpha1*exp(loglambda1)*td^(alpha1-1)*exp(-(exp(loglambda1)*td^alpha1 + exp(loglambda2)*td^alpha2)))
  lftr <- log(alpha2*exp(loglambda2)*tr^(alpha2-1)*exp(-(exp(loglambda1)*tr^alpha1 + exp(loglambda2)*tr^alpha2)))
  
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(lftd) + sum(lftr)
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}



#dat <- cfr_simdata
#params <- c(rep(0.001, 3), rep(0.01, 3), .05)

#dat <- cfr_simdata
#params <- c(rep(0.001, 3), rep(0.01, 3), .05)
LLerlang.ind <- function(params, data, compartments = 3){
  #parameters
  loglambda1 <- params[1:compartments]
  loglambda2 <- params[(compartments+1):(compartments+compartments)]
  alpha <- params[2*compartments+1]
  
  lambda1 <- exp(loglambda1)
  lambda2 <- exp(loglambda2)
  # data
  dat <- data
  c <- dat$c
  d <- dat$d
  r <- dat$r
  T <- dat$T
  L <- dat$L
  td <- dat$td
  tr <- dat$tr
  
  
  # individual level data
  # initialize expected values
  indiv_expected_risk = matrix(0, ncol=compartments, nrow=T)  # initialize expected_riskset to all 0s;
  indiv_expect_death <- vector(length = T)
  indiv_expect_recov<- vector(length = T)
  
  indiv_expected_risk[1,] <- 1
  indiv_expect_death[1] = indiv_expected_risk[1,1]*lambda1[1]
  indiv_expect_recov[1] = indiv_expected_risk[1,1]*lambda2[1]
  
  for (t in 2:T) {
    indiv_expected_risk[t,1] =  
      indiv_expected_risk[t-1,1]  - 
      indiv_expected_risk[t-1,1] * lambda1[1] -
      indiv_expected_risk[t-1,1] * lambda2[1]
    if (compartments>1) {
      for (i in 2:compartments) {
        indiv_expected_risk[t,i-1] =  indiv_expected_risk[t-1,i-1] - alpha *indiv_expected_risk[t-1,i-1];
        
        indiv_expected_risk[t,i] =  
          indiv_expected_risk[t-1,i] +
          alpha *indiv_expected_risk[t-1,i-1] -
          indiv_expected_risk[t-1,i] * lambda1[i] -
          indiv_expected_risk[t-1,i] * lambda2[i];
      }
    }
    
    #Accumulate deaths and recovereds                      
    indiv_expect_death[t] = indiv_expected_risk[t,1] * lambda1[1];
    indiv_expect_recov[t] = indiv_expected_risk[t,1] * lambda2[1];
    
    if(compartments>1) {
      for (i in 2:compartments) {
        indiv_expect_death[t] = indiv_expect_death[t-1] + indiv_expected_risk[t,i] * lambda1[i];
        indiv_expect_recov[t] = indiv_expect_recov[t-1] + indiv_expected_risk[t,i] * lambda2[i];
      }
    }
  }
  
  #aggregate data
  expected_deaths <- matrix(nrow = T, ncol = L)
  expected_recovereds <- matrix(nrow = T, ncol = L)
  for (j in 1:L) {
    expected_riskset <- matrix(nrow = T, ncol = compartments)
    
    expected_riskset[1,1] = c[1,j]+0.001
    expected_deaths[1,j] = expected_riskset[1,1] * lambda1[1]  + 0.0001 
    expected_recovereds[1,j] = expected_riskset[1,1] * lambda2[1]  + 0.0001
    
    #Initialize the riskset i ther compartments to 0 if they exist.
    if(compartments>1) {
      for (i in 2:compartments) {
        expected_riskset[1,i] = 0;
      }
    }
    
    for (t in 2:T) {
      #Move deaths and recoveries from last timestep out of first compartment.
      expected_riskset[t,1] =  expected_riskset[t-1,1] + c[t,j] - 
        expected_riskset[t-1,1] * lambda1[1]  - 
        expected_riskset[t-1,1] * lambda2[1] 
      
      if(compartments>1) {
        for (i in 2:compartments) {
          #remove the alpha from the expected risk set i-1
          expected_riskset[t,i-1] = expected_riskset[t-1,i-1] - alpha *expected_riskset[t-1,i-1]
          
          expected_riskset[t,i] =  expected_riskset[t-1,i] +
            alpha *expected_riskset[t-1,i-1] -
            expected_riskset[t-1,i] * lambda1[i] -
            expected_riskset[t-1,i] * lambda2[i] 
        }
      }                            
      
      #Accumulate deaths and recovereds                      
      expected_deaths[t,j] = expected_riskset[t,1] * lambda1[1]
      expected_recovereds[t,j] = expected_riskset[t,1] *  lambda2[1] 
      
      if(compartments>1) {
        for (i in 2:compartments) {
          expected_deaths[t,j] = expected_deaths[t-1,j] + expected_riskset[t,i] * lambda1[i] 
          expected_recovereds[t,j] =   expected_recovereds[t-1,j] + expected_riskset[t,i] * lambda2[i] 
        }
      }
      
    }
  }
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(log(indiv_expect_death[td])) + sum(log(indiv_expect_recov))
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}




