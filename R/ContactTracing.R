## COVID19 Contact Tracing

## Function to randomly draw a delay from 
## symptom onset to isolation, using a dist
## indexed from time of symptom onset
pullDelayOnset <- function(delay_par, n){
  return(rlnorm(n, meanlog=delay_par[1], sdlog=delay_par[2]))
}

## Function to randomly draw a delay from 
## symptom onset to isolation, using a dist
## indexed from time of infection
pullDelayInfection <- function(delay_par, inc_par, n){
  t_inf_iso <- rlnorm(n, meanlog=delay_par[1], sdlog=delay_par[2])
  t_inf_sym <- rlnorm(n, meanlog=inc_par[1], sdlog=inc_par[2])
  return(t = t_inf_iso - t_inf_sym)
}

## Function to randomly draw a delay from 
## symptom onset to isolation, using a dist
## indexed from time of isolation of prior generation
## (that is, delay = time from I0 isolation to I1 isolation)
pullDelayIsolation <- function(delay_par_G0, delay_par_G1, inf_par, inc_par, n){
  t_onset0_onset1 <- rlnorm(n, meanlog=inc_par[1], sdlog=inc_par[2]) + (rgamma(n, shape=inf_par[1], rate=inf_par[2]) - inf_par[3])
  t_onset0_iso0 <- rlnorm(n, meanlog=delay_par_G0[1], sdlog=delay_par_G0[2])
  t_iso0_iso1 <- rlnorm(n, meanlog=delay_par_G0[1], sdlog=delay_par_G0[2])
  t_onset1_iso1 <- t_iso0_iso1 + (t_onset0_iso0 - t_onset0_onset1)
  return(t_onset1_iso1)
}

## Function to determine proportion of infectious period
## reduced given distribution of infectiousness and a
## delay (time from symptom onset to detection)
## The distribution of infectiousness is from time of
## symptom onset, offset by inf_par[3]
calcInfProp <- function(delay, inf_par){
  return(pgamma(delay + inf_par[3], shape=inf_par[1], rate=inf_par[2]))
}

##' @name calcRho1
##' @description Function to calculate proportion of cases that can be detected in passive surveillance
##'
##' @param I0 the number of infected individuals in a generation 
##' @param p list of parameters with named parameters as elements
##'
##' @return rho1, probability of a given case being detected
##'
calcRho1 <- function(p, I0){
  if(!is.null(p[['testing_limit']])){
    if(is.null(p[['max_rho1']])){warning("assuming all cases potentially identifiable in passive surveillance")}else{
      if(p[['max_rho1']]>1){stop("Please specify max_rho1 as a probability in [0,1]")}}
    
    # note that this is a crude approximation and assumes that an
    # equal number of individuals would become infected per day
    rho1 <- min(max_rho1, p[['testing_limit']]/I0)
  }else{
    if(is.null(p[['rho1']])) stop("Must specify either testing_limit or rho1")
    rho1 <- p[['rho1']]
  }
  return(rho1)
}

##' @name calcRhoC
##' @description Function to calculate proportion of non-household cases that can be detected in passive surveillance
##'
##' @param C1_P_C the number of non-household contacts to follow in a generation
##' @param p list of parameters with named parameters as elements
##'
##' @return rhoC, probability of a non-household contact being successfully traced
##'
calcRhoC <- function(p, C1_H_C, C1_H_H){
  if(!is.null(p[['contact_limit']])){
    if(is.null(p[['max_rhoC']])){warning("assuming all non-household cases are potentially traceable")}else{
      if(p[['max_rhoC']]>1){stop("Please specify max_rhoC as a probability in [0,1]")}}
    

    if(is.null(p[['incl_hh']])) stop("Must specify whether household contacts are included in contact tracing limits")
    if(p[['incl_hh']]){
      rhoC <- min(max_rhoC, p[['contact_limit']]/(C1_P_C + CI_H_C))
    }else{
      rhoC <- min(max_rhoC, p[['contact_limit']]/C1_P_C)
    }
  }else{
    if(is.null(p[['rhoC']])) stop("Must specify either contact_limit or rhoC")
    rhoC <- p[['rhoC']]
  }
  return(rhoC)
}

## Function to perform simulation for given parameter set
##' @name run_tracing_sim
##' @description Function to run stochastic simulations of the contact tracing process from I0 - I2
##'
##' @param nsims number of simulations to run
##' @param I0 the number of initially infected individuals to begin each simulation
##' @param p list of parameters with named parameters as elements
##'
##' @return dataframe of final number infected in I0, I1, I2 for each simulation
##'
##'
run_tracing_sim <- function(nsims=1, I0, p){
  
  # Decide whether delays to isolation among household contacts are indexed from symptom onset, infection/exposure, or isolation of the prior generation
  if(p[['ind_H']]=="onset"){
    delayFuncHH <- pullDelayOnset
  }
  if(p[['ind_H']]=="infection"){
    delayFuncHH <- function(delay_par, inc_par = p[['inc_par']], n){return(pullDelayInfection(delay_par, inc_par, n))}  
  }
  if(p[['ind_H']]=="isolation"){
    delayFuncHH <- function(delay_par, delay_par_G0 = p[['delay_parP']], inc_par = p[['inc_par']], inf_par = p[['inf_par']], n){
      return(pullDelayIsolation(delay_par_G0, delay_par_G1=delay_par, inf_par, inc_par, n))}  
  }

  if(p[['ind_C']]=="onset"){
    delayFuncC <- pullDelayOnset
  }
  if(p[['ind_C']]=="infection"){
    delayFuncC <- function(delay_par, inc_par = p[['inc_par']], n){return(pullDelayInfection(delay_par, inc_par, n))}  
  }
  if(p[['ind_C']]=="isolation"){
    delayFuncC <- function(delay_par, delay_par_G0 = p[['delay_parP']], inc_par = p[['inc_par']], inf_par = p[['inf_par']], n){
      return(pullDelayIsolation(delay_par_G0, delay_par_G1=delay_par, inf_par, inc_par, n))}  
  }
  
  I0_vec <- rep(I0, nsims)

  rho1 <- calcRho1(p, I0_vec)
  
  # I0 split into I0_passive = I0*rho1, I0_undetected = I0*(1-rho1)
  I0_P <- rbinom(nsims, I0_vec, rho1)
  
    # I0_passive has contacts [C1_P] = I0*rho1*N
    ## this probably needs some sort of negative binom draw to match with R0 dist?
    C1_P <- rnbinom(length(I0_P), mu = p[['N']], size = p[['thetaN']])
    
      # [C1_HH] I0's HH contacts = [I0*rho1*N]*alphaH
      C1_P_HH <- rbinom(nsims, C1_P, p[['alphaH']])
      C1_P_C <- C1_P - C1_P_HH
      
      # Calculate detection rates for household + non-household contacts
      # Currently rhoH must be directly defined
      rhoC <- calcRhoC(p, C1_H_C, C1_P_HH)
      rhoH <- p[['rhoH']]
      
      # [C1_HH_T] total I0 HH contacts traced = [I0*rho1*N*alphaH]*rhoH
      C1_HH_T <- rbinom(nsims, C1_P_HH, rhoH)
      
      # [C1_C_T] total non I0 HH contacts traced = I0*rho1*N*(1-alphaH)*rhoC
      C1_C_T <- rbinom(nsims, C1_P_C, rhoC)
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*alphaH      
      ## this is currently assuming equal AR in I0's HH vs non HH contacts
      ## might also be better to convert to per-contact prob of infection and base on C1?
      I0_P_delay <- lapply(I0_P, function(x){pullDelayOnset(delay_par = p[['delay_parP']], n=x)})
      I0_P_gamma <- lapply(I0_P_delay, function(x){calcInfProp(x, p[['inf_par']])})
      I1_P <- sapply(1:length(I0_P), function(x) sum(rnbinom(I0_P[x], mu=I0_P_gamma[[x]]*p[['R0']], size=p[['theta']])))
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*(alphaH*nu / (alphaH*nu - alphaH + 1)
      I1_HH <- rbinom(nsims, I1_P, (p[['alphaH']]*p[['nu']])/(p[['alphaH']]*p[['nu']] - p[['alphaH']] +1))
      I1_HH_T <- rbinom(nsims, I1_HH, p[['rhoH']])
      I1_HH_U <- I1_HH - I1_HH_T
      
      # [I1] I0's non HH contacts that become infected = I0*rho1*gammaP*R0*(1-alphaH)      
      I1_C <- I1_P - I1_HH
      I1_C_T <- rbinom(nsims, I1_C, rhoC)
      I1_C_U <- I1_C - I1_C_T
      
          # [I2] infections from HH I1 untraced - [I0*rho1*gammaP*R0*alphaH]*(1-rhoH) * R0
          I2_HH_untraced <- sapply(I1_HH_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
          
          # [I2] infections from HH I1 traced - [I0*rho1*gammaP*R0*alphaH] * rhoH * gammaH*R0
          # to decide if overdispersion remains the same with scaled R0
          I1_HH_T_delay <- lapply(I1_HH_T, function(x){delayFuncHH(delay_par = p[['delay_parH']], n=x)})
          I1_HH_T_gamma <- lapply(I1_HH_T_delay, function(x){calcInfProp(x, p[['inf_par']])})
          I2_HH_traced <- sapply(1:length(I1_HH_T), function(x) sum(rnbinom(I1_HH_T[x], mu=I1_HH_T_gamma[[x]]*p[['R0']], size=p[['theta']])))
          
          # [I2] infections from non HH I1 untraced - [I0*rho1*gammaP*R0*(1-alphaH)] * (1-rhoC) * R0
          I2_C_untraced <- sapply(I1_C_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
      
          # [I2] infections from non HH I1 traced - [I0*rho1*gammaP*R0*(1-alphaH)] * rhoC * gammaC*R0
          # to decide if overdispersion remains the same with scaled R0
          I1_C_T_delay <- lapply(I1_C_T, function(x){delayFuncC(delay_par = p[['delay_parC']], n=x)})
          I1_C_T_gamma <- lapply(I1_C_T_delay, function(x){calcInfProp(x, p[['inf_par']])})
          I2_C_traced <- sapply(1:length(I1_C_T), function(x) sum(rnbinom(I1_C_T[x], mu=I1_C_T_gamma[[x]]*p[['R0']], size=p[['theta']])))
          
  # [I1] infections from undetected I0 = I0*R0*(1-rho1)
  I0_undet <- I0_vec - I0_P
  I1_U <- sapply(I0_undet, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
  
    # [I1] detected = I0*R0*(1-rho1)*rho1
    I1_U_P <- rbinom(nsims, I1_U, rho1)
  
    # [I1] undetected = I0*R0*(1-rho1)^2
    I1_U_U <- I1_U - I1_U_P
    
      # [I2] infections from undetected I1 = I0*(R0^2)*(1-rho1)^2
      I2_U_U <- sapply(I1_U_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
    
      # [I2] infections from detected I1 = I0*R0*(1-rho1)*rho1 * gammaP*R0
      I1_U_P_delay <- lapply(I1_U_P, function(x){pullDelayOnset(delay_par = p[['delay_parP']], n=x)})
      I1_U_P_gamma <- lapply(I1_U_P_delay, function(x){calcInfProp(x, p[['inf_par']])})
      I2_U_P <- sapply(1:length(I1_U_P), function(x) sum(rnbinom(I1_U_P[x], mu=I1_U_P_gamma[[x]]*p[['R0']], size=p[['theta']])))
    
  ## Book-keeping
  I1_total <- I1_P + I1_U
  I1_passive <- I1_P
  I1_undetected <- I1_U
  I1_passive_untraced <- I1_HH_U + I1_C_U
  I2_total <- I2_HH_traced + I2_HH_untraced + I2_C_traced + I2_C_untraced + I2_U_U + I2_U_P
  I2_passive <- I2_HH_traced + I2_HH_untraced + I2_C_traced + I2_C_untraced
  I2_passive_untraced <-  I2_HH_untraced + I2_C_untraced
  I2_undetected <- I2_U_P + I2_U_U
  C1_total <- C1_P
  C1_traced <- C1_HH_T + C1_C_T
  
  # to decide what other values to output
  return(out = data.frame(I0_vec, 
                          I1_total,
                          I1_passive,
                          I1_passive_untraced,
                          I1_undetected,
                          I2_total,
                          I2_passive,
                          I2_passive_untraced,
                          I2_undetected,
                          C1_total,
                          C1_traced))
}


## Function to calculate Re as ratio of simulated gen1 + gen2 cases
##' @name calcRe_sim
##' @description wrapper function to calculate Reff + 95% CI between two given generations
##'
##' @param gen1 vector of total infectious (infected) individuals in generation 1
##' @param gen2 vector of total infectious (infected) individuals in generation 2, infected by gen1
##'
##' @return named vector of mean, median, 2.5% and 97.5% quantile of Reff
##'
calcRe_sim <- function(gen1, gen2){
  xx <- gen2[gen1>0] / gen1[gen1>0]
  return(c(mean = mean(xx),
           median = median(xx),
           lo = quantile(xx, 0.025),
           hi = quantile(xx, 0.975)))
}

##' @name calcRe_exact
##' @description Function to calculate Re for a given strategy using closed form solutions
##'
##' @param p list of parameters with named parameters as elements
##'
##' @return named vector of Re_I0, Re_I1, Re_I1_passive, Re_I1_undetected
##'
calcRe_exact <- function(p){
  
  if(p[['ind_H']]=="onset"){
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2)
  }
  if(p[['ind_H']]=="infection"){
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
  }
  if(p[['ind_H']]=="isolation"){
    stop("not yet supported")
  }
  
  if(p[['ind_C']]=="onset"){
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2)
  }
  if(p[['ind_C']]=="infection"){
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
    #delayC = exp(p[['delay_parC']][1]) - exp(p[['inc_par']][1])
  }
  if(p[['ind_C']]=="isolation"){
    stop("not yet supported")
  }
  
  delayP = exp(p[['delay_parP']][1] + p[['delay_parP']][2]^2/2)

  # calculate average infection proportions
  gammaP = calcInfProp(delayP, p[['inf_par']])
  gammaH = calcInfProp(delayH, p[['inf_par']])
  gammaC = calcInfProp(delayC, p[['inf_par']])
  
  R0  <- p[['R0']]
  rho1 <- p[['rho1']]
  rhoH <- p[['rhoH']] 
  rhoC <- p[['rhoC']]
  alphaH <- p[['alphaH']]
  nu <- p[['nu']]
  
  v1 <-rho1*gammaP - rho1 + 1
  v2 <- alphaH*nu - alphaH + 1
  Re_I0 <- R0*v1
  Re_I1_undetected <- R0*v1
  Re_I1_passive <- R0/v2 * ( (alphaH*nu)*(rhoH*gammaH - rhoH + 1) + (1-alphaH)*(rhoC*gammaC - rhoC + 1) )
  Re_I1 <- (R0/v1) * ( (rho1*gammaP)/v2 * ( (alphaH*nu)*(rhoH*gammaH - rhoH +1) + (1-alphaH)*(rhoC*gammaC - rhoC + 1)) + (1-rho1)*v1)
  
  return(c(Re_I0=Re_I0, 
           Re_I1_passive=Re_I1_passive, 
           Re_I1_undetected=Re_I1_undetected,
           Re_I1=Re_I1))
}



##' @name calcRe_exact_asym
##' @description Function to calculate Re for a given strategy using closed 
##' form solutions, accounting for difference in symptomatic and asymptomatic
##'
##' @param p_asym list of parameters with named parameters as elements
##'
##' @return named vector of Re_I0, Re_I1, Re_I1_passive, Re_I1_undetected
##'
calcRe_exact_asym <- function(p_asym){
  
  p <- p_asym
  
  if(p[['ind_H']]=="onset"){
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2)
  }
  if(p[['ind_H']]=="infection"){
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
  }
  if(p[['ind_H']]=="isolation"){
    stop("not yet supported")
  }
  
  if(p[['ind_C']]=="onset"){
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2)
  }
  if(p[['ind_C']]=="infection"){
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
    #delayC = exp(p[['delay_parC']][1]) - exp(p[['inc_par']][1])
  }
  if(p[['ind_C']]=="isolation"){
    stop("not yet supported")
  }
  
  delayP_asym = exp(p[['delay_parP_asym']][1] + p[['delay_parP_asym']][2]^2/2)
  delayP_sym = exp(p[['delay_parP_sym']][1] + p[['delay_parP_sym']][2]^2/2)
  
  # calculate average infection proportions
  gammaPA = calcInfProp(delayP_asym, p[['inf_par']])
  gammaPS = calcInfProp(delayP_sym, p[['inf_par']])
  gammaH = calcInfProp(delayH, p[['inf_par']])
  gammaC = calcInfProp(delayC, p[['inf_par']])
  
  R0  <- p[['R0']]
  rhoPA <- p[['rhoPA']]
  rhoPS <- p[['rhoPS']]
  rhoH <- p[['rhoH']] 
  rhoC <- p[['rhoC']]
  alphaH <- p[['alphaH']]
  nu <- p[['nu']]
  tau <- p[['tau']]
  
  d <- alphaH*nu - alphaH + 1
  x <- tau*(1-rhoPS) + (1-tau)*(1-rhoPA)
  y <- tau*rhoPS*gammaPS + (1-tau)*rhoPA*gammaPA
  z <- x + y
  Re_I0 <- R0*z
  Re_I1_undetected <- R0*z
  Re_I1_passive <- (R0/d) * ( (alphaH*nu)*(rhoH*gammaH + (1 - rhoH)*z) + (1-alphaH)*(rhoC*gammaC + (1-rhoC)*Z) )
  Re_I1 <- R0* (x + (y/(z*d))*((alphaH*nu)*(rhoH*gammaH + (1 - rhoH)*z) + (1-alphaH)*(rhoC*gammaC + (1-rhoC)*Z)))
  
  return(c(Re_I0=Re_I0, 
           Re_I1_passive=Re_I1_passive, 
           Re_I1_undetected=Re_I1_undetected,
           Re_I1=Re_I1))
}
