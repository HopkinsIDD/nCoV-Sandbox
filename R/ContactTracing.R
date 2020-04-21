## COVID19 Contact Tracing
## updated 19 April 2020

## Function to randomly draw a delay to identification


## Function to determine proportion of infectious period
## reduced given distribution of infectiousness and a
## delay from infection to detection


## Functions to calculate number of tests - assuming
## all contacts get tested


## Function to calculate number of tests - assuming
## only non-symptomatic contacts get tested
  # determine asymptomatic + presymptomatic proportion at time t from infection



## Function to perform simulation for given parameter set
run_tracing_sim <- function(nsims=1, I0, R0, theta, N,
                            p1, gammaP,
                            alphaH, rhoH, gammaH, 
                            rhoC, gammaC){
  I0_vec <- rep(I0, nsims)

  # I0 split into I0_passive = I0*p1, I0_undetected = I0*(1-p1)
  I0_P <- rbinom(nsims, I0_vec, p1)
  
    # I0_passive has contacts [C1_P] = I0*p1*N
    ## this probably needs some sort of negative binom draw to match with R0 dist?
    C1_P <- I0_P * N
    
      # [C1_HH] I0's HH contacts = [I0*p1*N]*alphaH
      C1_P_HH <- rbinom(nsims, C1_P, alphaH)
      C1_P_C <- C1_P - C1_P_HH
      
      # [C1_HH_T] total I0 HH contacts traced = [I0*p1*N*alphaH]*rhoH
      C1_HH_T <- rbinom(nsims, C1_P_HH, rhoH)
      
      # [C1_C_T] total non I0 HH contacts traced = I0*p1*N*(1-alphaH)*rhoC
      C1_C_T <- rbinom(nsims, C1_P_C, rhoC)
      
      # [I1] I0's HH contacts that become infected = I0*p1*gammaP*R0*alphaH      
      ## this is currently assuming equal AR in I0's HH vs non HH contacts
      ## might also be better to convert to per-contact prob of infection and base on C1?
      I1_P <- sapply(I0_P, function(x) sum(rnbinom(x, mu=(1-gammaP)*R0, size=theta)))
      
      # [I1] I0's HH contacts that become infected = I0*p1*gammaP*R0*alphaH
      I1_HH <- rbinom(nsims, I1_P, alphaH)
      I1_HH_T <- rbinom(nsims, I1_HH, rhoH)
      I1_HH_U <- I1_HH - I1_HH_T
      
      # [I1] I0's non HH contacts that become infected = I0*p1*gammaP*R0*(1-alphaH)      
      I1_C <- I1_P - I1_HH
      I1_C_T <- rbinom(nsims, I1_C, rhoC)
      I1_C_U <- I1_C - I1_C_T
      
          # [I2] infections from HH I1 untraced - [I0*p1*gammaP*R0*alphaH]*(1-rhoH) * R0
          I2_HH_untraced <- sapply(I1_HH_U, function(x) sum(rnbinom(x, mu=R0, size=theta)))
          
          # [I2] infections from HH I1 traced - [I0*p1*gammaP*R0*alphaH] * rhoH * gammaH*R0
          # to decide if overdispersion remains the same with scaled R0
          I2_HH_traced <- sapply(I1_HH_T, function(x) sum(rnbinom(x, mu=(1-gammaH)*R0, size=theta)))
      
          # [I2] infections from non HH I1 untraced - [I0*p1*gammaP*R0*(1-alphaH)] * (1-rhoC) * R0
          I2_C_untraced <- sapply(I1_C_U, function(x) sum(rnbinom(x, mu=R0, size=theta)))
      
          # [I2] infections from non HH I1 traced - [I0*p1*gammaP*R0*(1-alphaH)] * rhoC * gammaC*R0
          # to decide if overdispersion remains the same with scaled R0
          I2_C_traced <- sapply(I1_C_T, function(x) sum(rnbinom(x, mu=(1-gammaC)*R0, size=theta)))
          
  # [I1] infections from undetected I0 = I0*R0*(1-p1)
  I0_undet <- I0_vec - I0_P
  I1_U <- sapply(I0_undet, function(x) sum(rnbinom(x, mu=R0, size=theta)))
  
    # [I1] detected = I0*R0*(1-p1)*p1
    I1_U_P <- rbinom(nsims, I1_U, p1)
  
    # [I1] undetected = I0*R0*(1-p1)^2
    I1_U_U <- I1_U - I1_U_P
    
      # [I2] infections from undetected I1 = I0*(R0^2)*(1-p1)^2
      I2_U_U <- sapply(I1_U_U, function(x) sum(rnbinom(x, mu=R0, size=theta)))
    
      # [I2] infections from detected I1 = I0*R0*(1-p1)*p1 * gammaP*R0
      I2_U_P <- sapply(I1_U_P, function(x) sum(rnbinom(x, mu=(1-gammaP)*R0, size=theta)))
    
  ## Book-keeping
  I1_total <- I1_P + I1_U
  I2_total <- I2_HH_traced + I2_HH_untraced + I2_C_traced + I2_C_untraced + I2_U_U + I2_U_P
  I2_CT <- I2_HH_traced + I2_HH_untraced + I2_C_traced + I2_C_untraced
  C1_T <- C1_HH_T + C1_C_T
  
  # to decide what other values to output
  return(out = data.frame(I0_vec, 
                          I1_total, 
                          I2_total,
                          I2_CT,
                          C1_T))
}

run_tracing_sim(nsims=10, I0=100, R0=2.5, theta=0.2, N=20, p1=0.1, gammaP=0.3, alphaH=0.5, rhoH=0.9, gammaH=0.1, rhoC=0.5, gammaC=0.5)
# gammaP/h/c = proportion of infectious period that remains (I think we had defined as reduction? should be (1-gammaP)*R0?)

