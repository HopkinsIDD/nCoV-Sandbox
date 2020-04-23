## COVID19 Contact Tracing
## updated 22 April 2020


## Function to perform simulation for given parameter set
run_tracing_sim <- function(nsims=1, I0, p){
  I0_vec <- rep(I0, nsims)

  # I0 split into I0_passive = I0*rho1, I0_undetected = I0*(1-rho1)
  I0_P <- rbinom(nsims, I0_vec, p['rho1'])
  
    # I0_passive has contacts [C1_P] = I0*rho1*N
    ## this probably needs some sort of negative binom draw to match with R0 dist?
    C1_P <- I0_P * p['N']
    
      # [C1_HH] I0's HH contacts = [I0*rho1*N]*alphaH
      C1_P_HH <- rbinom(nsims, C1_P, p['alphaH'])
      C1_P_C <- C1_P - C1_P_HH
      
      # [C1_HH_T] total I0 HH contacts traced = [I0*rho1*N*alphaH]*rhoH
      C1_HH_T <- rbinom(nsims, C1_P_HH, p['rhoH'])
      
      # [C1_C_T] total non I0 HH contacts traced = I0*rho1*N*(1-alphaH)*rhoC
      C1_C_T <- rbinom(nsims, C1_P_C, p['rhoC'])
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*alphaH      
      ## this is currently assuming equal AR in I0's HH vs non HH contacts
      ## might also be better to convert to per-contact prob of infection and base on C1?
      I1_P <- sapply(I0_P, function(x) sum(rnbinom(x, mu=p['gammaP']*p['R0'], size=p['theta'])))
      
      # [I1] I0's HH contacts that become infected = I0*rho1*gammaP*R0*alphaH
      I1_HH <- rbinom(nsims, I1_P, p['alphaH'])
      I1_HH_T <- rbinom(nsims, I1_HH, p['rhoH'])
      I1_HH_U <- I1_HH - I1_HH_T
      
      # [I1] I0's non HH contacts that become infected = I0*rho1*gammaP*R0*(1-alphaH)      
      I1_C <- I1_P - I1_HH
      I1_C_T <- rbinom(nsims, I1_C, p['rhoC'])
      I1_C_U <- I1_C - I1_C_T
      
          # [I2] infections from HH I1 untraced - [I0*rho1*gammaP*R0*alphaH]*(1-rhoH) * R0
          I2_HH_untraced <- sapply(I1_HH_U, function(x) sum(rnbinom(x, mu=p['R0'], size=p['theta'])))
          
          # [I2] infections from HH I1 traced - [I0*rho1*gammaP*R0*alphaH] * rhoH * gammaH*R0
          # to decide if overdispersion remains the same with scaled R0
          I2_HH_traced <- sapply(I1_HH_T, function(x) sum(rnbinom(x, mu=p['gammaH']*p['R0'], size=p['theta'])))
      
          # [I2] infections from non HH I1 untraced - [I0*rho1*gammaP*R0*(1-alphaH)] * (1-rhoC) * R0
          I2_C_untraced <- sapply(I1_C_U, function(x) sum(rnbinom(x, mu=p['R0'], size=p['theta'])))
      
          # [I2] infections from non HH I1 traced - [I0*rho1*gammaP*R0*(1-alphaH)] * rhoC * gammaC*R0
          # to decide if overdispersion remains the same with scaled R0
          I2_C_traced <- sapply(I1_C_T, function(x) sum(rnbinom(x, mu=p['gammaC']*p['R0'], size=p['theta'])))
          
  # [I1] infections from undetected I0 = I0*R0*(1-rho1)
  I0_undet <- I0_vec - I0_P
  I1_U <- sapply(I0_undet, function(x) sum(rnbinom(x, mu=p['R0'], size=p['theta'])))
  
    # [I1] detected = I0*R0*(1-rho1)*rho1
    I1_U_P <- rbinom(nsims, I1_U, p['rho1'])
  
    # [I1] undetected = I0*R0*(1-rho1)^2
    I1_U_U <- I1_U - I1_U_P
    
      # [I2] infections from undetected I1 = I0*(R0^2)*(1-rho1)^2
      I2_U_U <- sapply(I1_U_U, function(x) sum(rnbinom(x, mu=p['R0'], size=p['theta'])))
    
      # [I2] infections from detected I1 = I0*R0*(1-rho1)*rho1 * gammaP*R0
      I2_U_P <- sapply(I1_U_P, function(x) sum(rnbinom(x, mu=p['gammaP']*p['R0'], size=p['theta'])))
    
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
  
  Re_I0 <- p['R0']*((p['rho1'] - 1)*p['gammaP'] +1)
  Re_I1_passive <- p['R0'] *(p['rhoH']*p['alphaH']*(p['gammaH']-1) + p['rhoC']*(1-p['alphaH'])*(p['gammaC']-1)+1)
  Re_I1_undetected <- p['R0']*(p['rho1']*(p['gammaP']-1)+1)
  Re_I1 <- p['R0'] * ( (p['rho1']*p['gammaP']*(p['rhoH']*p['alphaH']*(p['gammaH']-1) + p['rhoC']*(1-p['alphaH'])*(p['gammaC']-1) +1)) /(p['rho1']*(p['gammaP']-1)+1) - p['rho1'] + 1)
    
  return(c(Re_I0=Re_I0, 
           Re_I1_passive=Re_I1_passive, 
           Re_I1_undetected=Re_I1_undetected,
           Re_I1=Re_I1))
}


### tests
nsims = 1000
I0 = 100
params <- c(R0 = 2.5,
            theta = 0.1,
            N = 10,
            rho1 = 0.1,
            gammaP = 0.3,
            alphaH = 0.5,
            rhoH = 0.9,
            gammaH = 0.1,
            rhoC = 0.5,
            gammaC = 0.5)

dat <- run_tracing_sim(nsims, I0, params)
calcRe_sim(gen1 = dat$I1_total, gen2 = dat$I2_total)
calcRe_sim(gen1 = dat$I1_passive, gen2 = dat$I2_passive)
calcRe_exact(params)


## TO DO
# delay function
# number tested

## Function to randomly draw a delay to identification


## Function to determine proportion of infectious period
## reduced given distribution of infectiousness and a
## delay from infection to detection


## Functions to calculate number of tests - assuming
## all contacts get tested


## Function to calculate number of tests - assuming
## only non-symptomatic contacts get tested
# determine asymptomatic + presymptomatic proportion at time t from infection

