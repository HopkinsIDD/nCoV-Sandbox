### Likelihood functions for MLE approach ###

##' likelihood function for exponential dist with only aggregate data
##'
##' @param params startinv values for parameters
##' @param data dataset including confirmed cases, deaths, recovereds on each day, T (# time points), L (# areas), td (vector of individual death times), and tr (vector of individual recovery times)
##'        in epi.curve with missing observations.
##' @return -1 * likelihood
#params <- init_params
LLexp.agg <- function(params, data, recovrate = "perfect"){
  # initialization parameters
  h1 <- exp(params[1])
  h2 <- exp(params[2])

  # format input data
 # dat <- data
  c <- data$c
  d <- data$d
  r <- data$r
  T <- data$T
  L <- data$L
  
  # initialize expected values
  expected_riskset = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0's
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_deaths to all 0's
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_recovereds to all 0's

  if(recovrate == "perfect"){
    # compute expected values
    for (j in 1:L) {  # Loops through countries (L)
      expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,j] * h1 + 0.001
      expected_recovereds[1,j] = expected_riskset[1,j] * h2 + 0.001
      for (t in 2:T) {  # Loops through times within each country
        expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                   (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
        expected_deaths[t,j] = expected_riskset[t,j] *  h1  +0.001
        expected_recovereds[t,j] = expected_riskset[t,j] *  h2 +0.001
      } 
    }
  }
  
  if(recovrate=="estimate"){
   #recover deteect rate
    dr <- expit(params[3])
    #is this the place with suboptimal recovery detection?
    w <- data$w
    for (j in 1:L) {  # Loops through countries (L)
      expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,j] * h1 + 0.001
      expected_recovereds[1,j] = expected_riskset[1,j] * h2 * dr^w[j] + 0.001
      for (t in 2:T) {  # Loops through times within each country
        expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                   (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
        expected_deaths[t,j] = expected_riskset[t,j] *  h1  +0.001
        expected_recovereds[t,j] = expected_riskset[t,j] *  h2 * dr^w[j] +0.001
      } 
    }
  }
  # compute expected values
  
  
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
LLexp.ind <- function(params, data, recovrate="perfect"){
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
  td <- dat$td  # individual-level data
  tr <- dat$tr  # individual-level data
  
  # initialize expected values
  expected_riskset = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_deaths = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  expected_recovereds = matrix(0, nrow=T, ncol=L)  # initialize expected_riskset to all 0s;
  
  #compute expected values for aggregate data
  
  if(recovrate=="perfect"){
    for (j in 1:L) {  # Loops through countries
      expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,j] * h1 +0.001
      expected_recovereds[1,j] = expected_riskset[1,j] * h2  +0.001
      for (t in 2:T) {  # Loops through all time points
        expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                   (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
        expected_deaths[t,j] = expected_riskset[t,j] *  h1  +0.001
        expected_recovereds[t,j] = expected_riskset[t,j] *  h2 +0.001
      } 
    }
  }
 
  
  if(recovrate=="estimate"){
    #recover deteect rate
    dr <- expit(params[3])
    #is this the place with suboptimal recovery detection?
    w <- data$w
    for (j in 1:L) {  # Loops through countries (L)
      expected_riskset[1,j] = c[1,j]  # risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,j] * h1 + 0.001
      expected_recovereds[1,j] = expected_riskset[1,j] * h2 * dr^w[j] + 0.001
      for (t in 2:T) {  # Loops through times within each country
        expected_riskset[t,j] = (expected_riskset[t-1,j] + c[t,j] - 
                                   (expected_riskset[t-1,j]*h1 +  (expected_riskset[t-1,j]* h2)))
        expected_deaths[t,j] = expected_riskset[t,j] *  h1  +0.001
        expected_recovereds[t,j] = expected_riskset[t,j] *  h2 * dr^w[j] +0.001
      } 
    }
  }
  # compute likelihood contributions for individual data: f(t)/F(tau)
  # lftd <- log(h1*exp(-(h1+h2)*td)) #f(t)
  # lftr <- log(h2*exp(-(h1+h2)*tr))
  
  lftd <- log(h1*exp(-(h1+h2)*td)/(h1/(h1+h2))) #f(t)/F(tau)
  lftr <- log(h2*exp(-(h1+h2)*tr)/(h2/(h1+h2)))
  
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


LLwbl.agg <- function(params, data, recovrate = "perfect"){
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
  
  if(recovrate=="perfect"){
    for (j in 1:L) {
      expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
      expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,1] * h1[1] +0.001
      expected_recovereds[1,j] = expected_riskset[1,1] * h2[1] +0.001
      
      for (t in 2:T) {
        expected_riskset[1,t] = c[t,j] 
        expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] +0.001
        expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1] +0.001
        
        for(k in 2:t){
          
          expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
          expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
          expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k]);
          
        }
      }
    }
  }
  
  if(recovrate=="estimate"){
    dr <- expit(params[5])
    w <- data$w
    for (j in 1:L) {
      expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
      expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,1] * h1[1] +0.001
      expected_recovereds[1,j] = expected_riskset[1,1] * h2[1] * dr^w[j] +0.001
      
      for (t in 2:T) {
        expected_riskset[1,t] = c[t,j] 
        expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] +0.001
        expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1] * dr^w[j] +0.001
        
        for(k in 2:t){
          
          expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
          expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
          expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k] * dr^w[j]);
          
        }
      }
    }
  }
  
  

  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds)))
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}


LLwbl.ind <- function(params, data, recovrate = "perfect"){
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
  
  if(recovrate=="perfect"){
    for (j in 1:L) {
      expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
      expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,1] * h1[1] +0.001
      expected_recovereds[1,j] = expected_riskset[1,1] * h2[1] +0.001
      
      for (t in 2:T) {
        expected_riskset[1,t] = c[t,j] 
        expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] +0.001
        expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1] +0.001
        
        for(k in 2:t){
          
          expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
          expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
          expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k]);
          
        }
      }
    }
  }
  
  if(recovrate=="estimate"){
    dr <- expit(params[5])
    w <- data$w
    for (j in 1:L) {
      expected_riskset = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
      expected_riskset[1,1] = c[1,j] #risk set on day 1 at inf duration 1 is just c1
      expected_deaths[1,j] = expected_riskset[1,1] * h1[1] +0.001
      expected_recovereds[1,j] = expected_riskset[1,1] * h2[1] * dr^w[j] +0.001
      
      for (t in 2:T) {
        expected_riskset[1,t] = c[t,j] 
        expected_deaths[t,j] = expected_riskset[1,t] *  h1[1] +0.001
        expected_recovereds[t,j] = expected_riskset[1,t] *  h2[1] * dr^w[j] +0.001
        
        for(k in 2:t){
          
          expected_riskset[k,t] = (expected_riskset[k-1,t-1] - (expected_riskset[k-1,t-1]*h1[k-1] +  (expected_riskset[k-1,t-1]* h2[k-1])));
          expected_deaths[t,j] =  expected_deaths[t,j] + (expected_riskset[k,t] * h1[k]);
          expected_recovereds[t,j] = expected_recovereds[t,j] + (expected_riskset[k,t] * h2[k] * dr^w[j]);
          
        }
      }
    }
  }
  
  #compute individual leve contributions
  
  indiv_expected_risk = matrix(0, nrow=T, ncol=T)  # initialize expected_riskset to all 0s;
  indiv_expect_death = rep(0, T)
  indiv_expect_recov = rep(0, T)
  indiv_expected_risk[1,1] = 1 #risk set on day 1 at inf duration 1 is just c1
  indiv_expect_death[1] = indiv_expected_risk[1,1] * h1[1] +0.001
  indiv_expect_recov[1] = indiv_expected_risk[1,1] * h2[1] +0.001
  
  for (t in 2:T) {
    indiv_expected_risk[1,t] = c[t,j] 
    indiv_expect_death[t] = indiv_expected_risk[1,t] *  h1[1]  +0.001
    indiv_expect_recov[t] = indiv_expected_risk[1,t] *  h2[1]+0.001
    
    for(k in 2:t){
      
      indiv_expected_risk[k,t] = (indiv_expected_risk[k-1,t-1] - (indiv_expected_risk[k-1,t-1]*h1[k-1] +  (indiv_expected_risk[k-1,t-1]* h2[k-1])));
      indiv_expect_death[t] =  indiv_expect_death[t] + (indiv_expected_risk[k,t] * h1[k]);
      indiv_expect_recov[t] = indiv_expect_recov[t] + (indiv_expected_risk[k,t] * h2[k]);
      
    }
  }
  
  #compute likelihood contributionss for individual data
  #f(t)
  # lftd <- log(alpha1*exp(loglambda1)*td^(alpha1-1)*exp(-(exp(loglambda1)*td^alpha1 + exp(loglambda2)*td^alpha2)))
  # lftr <- log(alpha2*exp(loglambda2)*tr^(alpha2-1)*exp(-(exp(loglambda1)*tr^alpha1 + exp(loglambda2)*tr^alpha2)))
  
  #f(t)/F(tau)
  lftd <- log(alpha1*exp(loglambda1)*td^(alpha1-1)*exp(-(exp(loglambda1)*td^alpha1 + exp(loglambda2)*td^alpha2)) / sum(indiv_expect_death))
  lftr <- log(alpha2*exp(loglambda2)*tr^(alpha2-1)*exp(-(exp(loglambda1)*tr^alpha1 + exp(loglambda2)*tr^alpha2)) / sum(indiv_expect_recov))
  
  
  
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(lftd) + sum(lftr)
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}




#data <- cfr_data
#params <- c(rep(log(0.001), 3), rep(log(0.01), 3), log(.05))
LLerlang.ind <- function(params, data, compartments = 3, recovrate="perfect"){
  #parameters
  loglambda1 <- params[1:compartments]
  loglambda2 <- params[(compartments+1):(compartments+compartments)]
  alpha <- exp(params[2*compartments+1])
 
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
#indiv_expect_death <- vector(length = T)
#indiv_expect_recov<- vector(length = T)
indiv_expect_death <-rep(0, T)
indiv_expect_recov<- rep(0, T)

indiv_expected_risk[1,1] <- 1
indiv_expect_death[1] = indiv_expected_risk[1,1]*lambda1[1]
indiv_expect_recov[1] = indiv_expected_risk[1,1]*lambda2[1]

for (t in 2:T) {
  indiv_expected_risk[t,1] =  
      indiv_expected_risk[t-1,1]  - 
      indiv_expected_risk[t-1,1] * lambda1[1] -
      indiv_expected_risk[t-1,1] * lambda2[1]

  
  if (compartments>1) {
    for (i in 2:compartments) {
      indiv_expected_risk[t,i-1] =  indiv_expected_risk[t,i-1] - alpha *indiv_expected_risk[t-1,i-1]
                                        
      
      indiv_expected_risk[t,i] =  
          indiv_expected_risk[t-1,i] +
          alpha *indiv_expected_risk[t-1,i-1] -
          indiv_expected_risk[t-1,i] * lambda1[i] -
          indiv_expected_risk[t-1,i] * lambda2[i]
  
      
    }
  }
  
  #Accumulate deaths and recovereds                      
  indiv_expect_death[t] = indiv_expected_risk[t,1] * lambda1[1] +0.001
  indiv_expect_recov[t] = indiv_expected_risk[t,1] * lambda2[1]+0.001
  
  if(compartments>1) {
    for (i in 2:compartments) {
      indiv_expect_death[t] = indiv_expect_death[t] + indiv_expected_risk[t,i] * lambda1[i];
      indiv_expect_recov[t] = indiv_expect_recov[t] + indiv_expected_risk[t,i] * lambda2[i];
    }
  }
}

#aggregate data
if(recovrate=="perfect"){
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
          #expected_riskset[t,i-1] = expected_riskset[t-1,i-1] - alpha *expected_riskset[t-1,i-1]
          expected_riskset[t,i-1] = expected_riskset[t,i-1] - alpha *expected_riskset[t-1,i-1]
          
          expected_riskset[t,i] =  expected_riskset[t-1,i] +
            alpha *expected_riskset[t-1,i-1] -
            expected_riskset[t-1,i] * lambda1[i] -
            expected_riskset[t-1,i] * lambda2[i] 
        }
      }                            
      
      #Accumulate deaths and recovereds                      
      expected_deaths[t,j] = expected_riskset[t,1] * lambda1[1] +0.001
      expected_recovereds[t,j] = expected_riskset[t,1] *  lambda2[1] +0.001
      
      if(compartments>1) {
        for (i in 2:compartments) {
          expected_deaths[t,j] = expected_deaths[t,j] + expected_riskset[t,i] * lambda1[i] 
          expected_recovereds[t,j] =   expected_recovereds[t,j] + expected_riskset[t,i] * lambda2[i] 
        }
      }
      
    }
  }
}



if(recovrate=="estimate"){
  dr <- expit(params[2*compartments+2])
  w <- data$w
  expected_deaths <- matrix(nrow = T, ncol = L)
  expected_recovereds <- matrix(nrow = T, ncol = L)
  for (j in 1:L) {
    expected_riskset <- matrix(nrow = T, ncol = compartments)
    
    expected_riskset[1,1] = c[1,j]+0.001
    expected_deaths[1,j] = expected_riskset[1,1] * lambda1[1]  + 0.0001 
    expected_recovereds[1,j] = expected_riskset[1,1] * lambda2[1]  * dr^w[j] + 0.0001
    
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
          #expected_riskset[t,i-1] = expected_riskset[t-1,i-1] - alpha *expected_riskset[t-1,i-1]
          expected_riskset[t,i-1] = expected_riskset[t,i-1] - alpha *expected_riskset[t-1,i-1]
          
          expected_riskset[t,i] =  expected_riskset[t-1,i] +
            alpha *expected_riskset[t-1,i-1] -
            expected_riskset[t-1,i] * lambda1[i] -
            expected_riskset[t-1,i] * lambda2[i] 
        }
      }                            
      
      #Accumulate deaths and recovereds                      
      expected_deaths[t,j] = expected_riskset[t,1] * lambda1[1] + 0.001
      expected_recovereds[t,j] = expected_riskset[t,1] *  lambda2[1]  * dr^w[j]  +0.001
      
      if(compartments>1) {
        for (i in 2:compartments) {
          expected_deaths[t,j] = expected_deaths[t,j] + expected_riskset[t,i] * lambda1[i] 
          expected_recovereds[t,j] =   expected_recovereds[t,j] + expected_riskset[t,i] * lambda2[i]  * dr^w[j] 
        }
      }
      
    }
  }
}

#likelihood
#LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(log(indiv_expect_death[td])) + sum(log(indiv_expect_recov[tr])) # indiv contribution = f(t)
LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) + sum(log(indiv_expect_death[td]/sum(indiv_expect_death))) + sum(log(indiv_expect_recov[tr]/sum(indiv_expect_recov))) #f(t)/F(tau)
#optim minimizes rather than maximizes, so take negative here
return(-LL)
}

LLerlang.agg <- function(params, data, compartments = 3, recovrate="perfect"){
  #parameters
  loglambda1 <- params[1:compartments]
  loglambda2 <- params[(compartments+1):(compartments+compartments)]
  alpha <- exp(params[2*compartments+1])
  
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
 
  
  #aggregate data
  if(recovrate=="perfect"){
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
            #expected_riskset[t,i-1] = expected_riskset[t-1,i-1] - alpha *expected_riskset[t-1,i-1]
            expected_riskset[t,i-1] = expected_riskset[t,i-1] - alpha *expected_riskset[t-1,i-1]
            
            expected_riskset[t,i] =  expected_riskset[t-1,i] +
              alpha *expected_riskset[t-1,i-1] -
              expected_riskset[t-1,i] * lambda1[i] -
              expected_riskset[t-1,i] * lambda2[i] 
          }
        }                            
        
        #Accumulate deaths and recovereds                      
        expected_deaths[t,j] = expected_riskset[t,1] * lambda1[1] +0.001
        expected_recovereds[t,j] = expected_riskset[t,1] *  lambda2[1] +0.001
        
        if(compartments>1) {
          for (i in 2:compartments) {
            expected_deaths[t,j] = expected_deaths[t,j] + expected_riskset[t,i] * lambda1[i] 
            expected_recovereds[t,j] =   expected_recovereds[t,j] + expected_riskset[t,i] * lambda2[i] 
          }
        }
        
      }
    }
  }
  
  
  
  if(recovrate=="estimate"){
    dr <- expit(params[2*compartments+2])
    w <- data$w
    expected_deaths <- matrix(nrow = T, ncol = L)
    expected_recovereds <- matrix(nrow = T, ncol = L)
    for (j in 1:L) {
      expected_riskset <- matrix(nrow = T, ncol = compartments)
      
      expected_riskset[1,1] = c[1,j]+0.001
      expected_deaths[1,j] = expected_riskset[1,1] * lambda1[1]  + 0.0001 
      expected_recovereds[1,j] = expected_riskset[1,1] * lambda2[1]  * dr^w[j] + 0.0001
      
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
            #expected_riskset[t,i-1] = expected_riskset[t-1,i-1] - alpha *expected_riskset[t-1,i-1]
            expected_riskset[t,i-1] = expected_riskset[t,i-1] - alpha *expected_riskset[t-1,i-1]
            
            expected_riskset[t,i] =  expected_riskset[t-1,i] +
              alpha *expected_riskset[t-1,i-1] -
              expected_riskset[t-1,i] * lambda1[i] -
              expected_riskset[t-1,i] * lambda2[i] 
          }
        }                            
        
        #Accumulate deaths and recovereds                      
        expected_deaths[t,j] = expected_riskset[t,1] * lambda1[1] + 0.001
        expected_recovereds[t,j] = expected_riskset[t,1] *  lambda2[1]  * dr^w[j]  +0.001
        
        if(compartments>1) {
          for (i in 2:compartments) {
            expected_deaths[t,j] = expected_deaths[t,j] + expected_riskset[t,i] * lambda1[i] 
            expected_recovereds[t,j] =   expected_recovereds[t,j] + expected_riskset[t,i] * lambda2[i]  * dr^w[j] 
          }
        }
        
      }
    }
  }
  
  #likelihood
  LL <- sum((-expected_deaths+d*log(expected_deaths)) + (-expected_recovereds + r*log(expected_recovereds))) #+ sum(log(indiv_expect_death[td])) + sum(log(indiv_expect_recov[tr]))
  
  #optim minimizes rather than maximizes, so take negative here
  return(-LL)
}



#to fix individual level data problem, condition on having died

# to do othis, divide f(t)/F(tau or max individual follow-up) --> check
