## COVID19 Contact Tracing
## updated 27 April 2020

## Function to randomly draw a delay from 
## symptom onset to isolation, using a dist
## indexed from time of symptom onset
pullDelay <- function(delay_par, n){
  return(rlnorm(n, meanlog=delay_par[1], sdlog=delay_par[2]))
}

## Function to randomly draw a delay from 
## symptom onset to isolation, using a dist
## indexed from time of infection
pullDelayInfection <- function(delay_par, inc_par, n){
  t_onset_iso <- rlnorm(n, meanlog=delay_par[1], sdlog=delay_par[2])
  t_onset_sym <- rlnorm(n, meanlog=inc_par[1], sdlog=inc_par[2])
  return(t = t_onset_iso - t_onset_sym)
}

## Function to determine proportion of infectious period
## reduced given distribution of infectiousness and a
## delay (time from symptom onset to detection)
## The distribution of infectiousness is from time of
## symptom onset, offset by inf_par[3]
calcInfProp <- function(delay, inf_par){
  return(pgamma(delay + inf_par[3], shape=inf_par[1], rate=inf_par[2]))
}

## Function to perform simulation for given parameter set
run_tracing_sim <- function(nsims=1, I0, p){
  
  # Decide whether delays to isolation are indexed from symptom onset or infection
  if(p[['ind_onset']]){delayFunc <- pullDelay}else{delayFunc <- function(delay_par, inc_par = p[['inc_par']], n){return(pullDelayInfection(delay_par, inc_par, n))}}
  
  I0_vec <- rep(I0, nsims)

  # I0 split into I0_passive = I0*rho1, I0_undetected = I0*(1-rho1)
  I0_P <- rbinom(nsims, I0_vec, p[['rho1']])
  
    # I0_passive has contacts [C1_P] = I0*rho1*N
    ## this probably needs some sort of negative binom draw to match with R0 dist?
    C1_P <- rnbinom(length(I0_P), mu = p[['N']], size = p[['thetaN']])
    
      # [C1_HH] I0's HH contacts = [I0*rho1*N]*alphaH
      C1_P_HH <- rbinom(nsims, C1_P, p[['alphaH']])
      C1_P_C <- C1_P - C1_P_HH
      
      # [C1_HH_T] total I0 HH contacts traced = [I0*rho1*N*alphaH]*rhoH
      C1_HH_T <- rbinom(nsims, C1_P_HH, p[['rhoH']])
      
      # [C1_C_T] total non I0 HH contacts traced = I0*rho1*N*(1-alphaH)*rhoC
      C1_C_T <- rbinom(nsims, C1_P_C, p[['rhoC']])
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*alphaH      
      ## this is currently assuming equal AR in I0's HH vs non HH contacts
      ## might also be better to convert to per-contact prob of infection and base on C1?
      I0_P_delay <- lapply(I0_P, function(x){delayFunc(delay_par = p[['delay_parP']], n=x)})
      I0_P_gamma <- lapply(I0_P_delay, function(x){calcInfProp(x, p[['inf_par']])})
      I1_P <- sapply(1:length(I0_P), function(x) sum(rnbinom(I0_P[x], mu=I0_P_gamma[[x]]*p[['R0']], size=p[['theta']])))
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*(alphaH*nu / (alphaH*nu - alphaH + 1)
      I1_HH <- rbinom(nsims, I1_P, (p[['alphaH']]*p[['nu']])/(p[['alphaH']]*p[['nu']] - p[['alphaH']] +1))
      I1_HH_T <- rbinom(nsims, I1_HH, p[['rhoH']])
      I1_HH_U <- I1_HH - I1_HH_T
      
      # [I1] I0's non HH contacts that become infected = I0*rho1*gammaP*R0*(1-alphaH)      
      I1_C <- I1_P - I1_HH
      I1_C_T <- rbinom(nsims, I1_C, p[['rhoC']])
      I1_C_U <- I1_C - I1_C_T
      
          # [I2] infections from HH I1 untraced - [I0*rho1*gammaP*R0*alphaH]*(1-rhoH) * R0
          I2_HH_untraced <- sapply(I1_HH_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
          
          # [I2] infections from HH I1 traced - [I0*rho1*gammaP*R0*alphaH] * rhoH * gammaH*R0
          # to decide if overdispersion remains the same with scaled R0
          I1_HH_T_delay <- lapply(I1_HH_T, function(x){delayFunc(delay_par = p[['delay_parH']], n=x)})
          I1_HH_T_gamma <- lapply(I1_HH_T_delay, function(x){calcInfProp(x, p[['inf_par']])})
          I2_HH_traced <- sapply(1:length(I1_HH_T), function(x) sum(rnbinom(I1_HH_T[x], mu=I1_HH_T_gamma[[x]]*p[['R0']], size=p[['theta']])))
          
          # [I2] infections from non HH I1 untraced - [I0*rho1*gammaP*R0*(1-alphaH)] * (1-rhoC) * R0
          I2_C_untraced <- sapply(I1_C_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
      
          # [I2] infections from non HH I1 traced - [I0*rho1*gammaP*R0*(1-alphaH)] * rhoC * gammaC*R0
          # to decide if overdispersion remains the same with scaled R0
          I1_C_T_delay <- lapply(I1_C_T, function(x){delayFunc(delay_par = p[['delay_parC']], n=x)})
          I1_C_T_gamma <- lapply(I1_C_T_delay, function(x){calcInfProp(x, p[['inf_par']])})
          I2_C_traced <- sapply(1:length(I1_C_T), function(x) sum(rnbinom(I1_C_T[x], mu=I1_C_T_gamma[[x]]*p[['R0']], size=p[['theta']])))
          
  # [I1] infections from undetected I0 = I0*R0*(1-rho1)
  I0_undet <- I0_vec - I0_P
  I1_U <- sapply(I0_undet, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
  
    # [I1] detected = I0*R0*(1-rho1)*rho1
    I1_U_P <- rbinom(nsims, I1_U, p[['rho1']])
  
    # [I1] undetected = I0*R0*(1-rho1)^2
    I1_U_U <- I1_U - I1_U_P
    
      # [I2] infections from undetected I1 = I0*(R0^2)*(1-rho1)^2
      I2_U_U <- sapply(I1_U_U, function(x) sum(rnbinom(x, mu=p[['R0']], size=p[['theta']])))
    
      # [I2] infections from detected I1 = I0*R0*(1-rho1)*rho1 * gammaP*R0
      I1_U_P_delay <- lapply(I1_U_P, function(x){delayFunc(delay_par = p[['delay_parP']], n=x)})
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
calcRe_sim <- function(gen1, gen2){
  xx <- gen2[gen1>0] / gen1[gen1>0]
  return(c(mean = mean(xx),
           median = median(xx),
           lo = quantile(xx, 0.025),
           hi = quantile(xx, 0.975)))
}


## function to estimate Re directly from parameters
calcRe_exact <- function(p){
  
  # calculate average delays
  if(p[['ind_onset']]){
    delayP = exp(p[['delay_parP']][1] + p[['delay_parP']][2]^2/2)
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2)
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2)
  }else{
    delayP = exp(p[['delay_parP']][1] + p[['delay_parP']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
    delayH = exp(p[['delay_parH']][1] + p[['delay_parH']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
    delayC = exp(p[['delay_parC']][1] + p[['delay_parC']][2]^2/2) - exp(p[['inc_par']][1] + p[['inc_par']][2]^2/2)
  }
  
  # calculate average infection proportions
  gammaP = calcInfProp(delayP, p[['inf_par']])
  gammaH = calcInfProp(delayH, p[['inf_par']])
  gammaC = calcInfProp(delayC, p[['inf_par']])
  
  Re_I0 <- p[['R0']]*((gammaP - 1)*p[['rho1']] +1)
  Re_I1_passive <- p[['R0']]*(p[['rhoH']]*p[['alphaH']]*(gammaH-1) + p[['rhoC']]*(1-p[['alphaH']])*(gammaC-1)+1)
  Re_I1_undetected <- p[['R0']]*(p[['rho1']]*(gammaP-1)+1)
  Re_I1 <- p[['R0']] * ( (p[['rho1']]*gammaP*(p[['rhoH']]*p[['alphaH']]*(gammaH-1) + p[['rhoC']]*(1-p[['alphaH']])*(gammaC-1) +1)) /(p[['rho1']]*(gammaP-1)+1) - p[['rho1']] + 1)
    
  return(c(Re_I0=Re_I0, 
           Re_I1_passive=Re_I1_passive, 
           Re_I1_undetected=Re_I1_undetected,
           Re_I1=Re_I1))
}


### tests
nsims = 10
I0 = 100
params <- list(R0 = 2.5,
               theta = 0.1,
               N = 10,
               thetaN = 0.1,
               alphaH = 0.5,
               rho1 = 0.1,
               rhoH = 0.9,
               rhoC = 0.5,
               nu = 1,
               inf_par = c(2.1157790, 0.6898583, 2.3066912),
               inc_par = c(1.621, 0.418),
               delay_parP = c(1.22, 0.78),
               delay_parH = c(0.77, 0.67),
               delay_parC = c(0.77, 0.67),
               ind_onset = TRUE)
 
dat <- run_tracing_sim(nsims, I0, params)
calcRe_sim(gen1 = dat$I1_total, gen2 = dat$I2_total)
calcRe_sim(gen1 = dat$I1_passive, gen2 = dat$I2_passive)
calcRe_exact(params)


## TO DO
# add in thresholds to calculate rho1/rhoH/rhoC
# third generation
# check math with nu
# number tested

## Functions to calculate number of tests - assuming
## all contacts get tested


## Function to calculate number of tests - assuming
## only non-symptomatic contacts get tested
# determine asymptomatic + presymptomatic proportion at time t from infection


