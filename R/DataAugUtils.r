


##'
##' Function to plot the estimated and reported case counts
##'
##' @param conf_cases Confirmed case data from JHU CSSE
##' @param incid_ests Estimated incidence 
##' @param locations  Locations to plot, can be a vector
##' 
##' @return a corrected version of the inferred incidence data corrected for reporting changes in Hubei.
##'
plot_incidence_ests_report_b <- function(conf_cases=jhucsse, incid_ests=incidence_data, locations="Hubei"){
     
     incid_ests <- incid_ests %>% mutate(Incidence = ceiling(Incidence)) %>% filter(Province_State %in% locations)
     incid_ests$Date <- as.Date(incid_ests$Date, "%m/%d/%Y")
     
     # Make conf_cases daily
     # Get daily calculated incidence (from reporting)
     conf_cases_daily <- conf_cases %>% filter(Province_State %in% locations & !is.na(Confirmed)) %>%
          mutate(Date = as.Date(Update)) %>% group_by(Province_State, Date) %>% filter(Update == max(Update, na.rm=TRUE)) %>% ungroup()
     
     conf_cases_daily <- conf_cases_daily %>% group_by(Province_State) %>% arrange(Date) %>% mutate(Incidence = diff(c(0, Confirmed))) %>% ungroup()
     
     # Merge in reported incidence
     incid_data_ <- left_join(incid_ests, conf_cases_daily %>% rename(Incid_rep = Incidence), by=c("Province_State", "Date")) %>% as_tibble()
     
     ggplot(incid_data_, aes(x=Date, y=Incidence)) + 
          geom_bar(stat = "identity", fill="maroon") +
          geom_point(aes(x=Date, y=Incid_rep), col="navyblue", size=3, alpha=0.5) +
          coord_cartesian(xlim=c(as.Date("2020-01-15"), max(as.Date(incid_data_$Date)))) +
          theme_classic() + theme(axis.text.x=element_text(size=11),
                                  axis.text.y=element_text(size=11),
                                  axis.title.x=element_text(size=13, margin = margin(t = 15)),
                                  axis.title.y=element_text(size=13, margin = margin(r = 15)),
                                  legend.position='none')
}



##' Function to simulate incubation period using log-normal distribution
##' 
##' @param n number of replicates
##' @param mu mean of lognormal dist
##' @param sigma standard devoiation of lognormal dist 
##' @param mu_95low lower 95% confidence interval around mu
##' @param mu_95high upper 95% confidence interval around mu
##' @param sigma_95low lower 95% confidence interval around sigma
##' @param sigma_95high upper 95% confidence interval around sigma
##' @param n_obs number of observations used to estimate confidence bounds of mu and sigma
##' 
##' @return scalar or vector of length n giving the random values of incubation period in days
##' 
sim_inf_period <- function(n, 
                           mu=1.621,
                           sigma=0.418,
                           mu_95low=1.504,
                           mu_95high=1.755,
                           sigma_95low=0.271,
                           sigma_95high=0.542,
                           n_obs=181
) {
     
     if (!is.null(c(mu_95low, mu_95high, sigma_95low, sigma_95high))) {
          
          # Get standard deviation from confidence intervals
          mu_sd <- log(sqrt(n_obs) * (exp(mu_95high) - exp(mu_95low)) / 3.92)
          sigma_sd <- log(sqrt(n_obs) * (exp(sigma_95high) - exp(sigma_95low)) / 3.92)
          
          pars <- cbind(truncnorm::rtruncnorm(n, a=mu_95low, b=mu_95high, mean=mu, sd=mu_sd),
                        truncnorm::rtruncnorm(n, a=sigma_95low, b=sigma_95high, mean=sigma, sd=sigma_sd))
          
          return(apply(pars, 1, function(x) rlnorm(1, x[1], x[2])))
          
     } else {
          
          return(rlnorm(n, mu, sigma))
     }
}


##' Downloads and cleans public linelist from Mortitz Kramer's googlesheets
##' 
##' @param save_csv save as csv (default=TRUE)
##' 
##' @return data frame 
##' 
read_MK_direct <- function(google=TRUE,
                           save_csv=TRUE
){
     
     require(googlesheets4)
     
     linelist <- 'https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/htmlview#'
     linelist_Hubei <- data.frame(read_sheet(linelist, sheet='Hubei'), stringsAsFactors = F)
     linelist_nonHubei <- data.frame(read_sheet(linelist, sheet='outside_Hubei'), stringsAsFactors = F)
     
     message('Combining and cleaning')
     d <- rbind(linelist_Hubei, 
                linelist_nonHubei[,names(linelist_nonHubei) != "data_moderator_initials"])
     
     # Keep only exact ages
     #d$age <- unlist(lapply(d$age, function(x) {
     #     tmp <- is.numeric(x)
     #     if (tmp == TRUE) {
     #          return(ceiling(x))
     #     } else {
     #          return(NA)
     #    }
     #}))
     
     d$sex[d$sex == 'N/A'] <- NA
     d <- d[!is.na(d$province) & d$province != 'China',] # observations missing province dropped

     d <- data.frame(age=as.numeric(d$age),
                      gender=factor(d$sex),
                      city=factor(d$city),
                      province=factor(d$province),
                      country=factor(d$country),
                      date_onset=as.Date(d$date_onset_symptoms, format='%d.%m.%Y'),
                      date_confirmation=as.Date(d$date_confirmation, format='%d.%m.%Y'),
                      date_admission=as.Date(d$date_admission_hospital, format='%d.%m.%Y'),
                      admitted=as.integer(!is.na(d$date_admission_hospital)),
                      death=as.integer(d$outcome == 'died'))
     
     d$time_to_care <- d$date_admission - d$date_onset
     d$time_to_care[d$time_to_care < 0] <- NA 
     d$time_to_conf <- d$date_confirmation - d$date_onset
     
     if(save_csv) {
          
          write_csv(d, paste0("./data/Linelist_MK_", format(Sys.Date(), "%Y_%m_%d"),".csv")) 
          message('Saved to .csv')
     } 
     
     return(d)
}

##' Fit Negative Binomial time to confirmation model
##' 
##' The function fits a Negative Binomial bayesian model to estimate the delay from time of symptom onset to time of case confirmation at the population-level (ignores spatial variation among provinces). 
##' The model uses a negative binomial error structure with probability p and dispersion parameter r. Model fit to data using JAGS and 'runjags' pacakge. Defualt MCMC params should be increased for full run.
##' 
##' @param data data frame of linelist data produced by the 'read_MK_direct' function
##' @param n_chains number of chains to run (default=4)
##' @param n_burn number burn in samples (default=1000)
##' @param n_samps number of samples (default=1000)
##' @param n_thin interval to thin samples (default=1)
##' 
##' @return mcmc model object
##' 

fit_negbin_time_to_conf <- function(data,
                                    n_chains=4, # number chains
                                    n_burn=1000, # burn in
                                    n_samp=1000, # samples
                                    n_thin=1 # thin
) {
     
     require(runjags)
     
     Y <- data$time_to_conf[!is.na(data$time_to_conf)]
     
     jags_data <- list(
          Y=Y, # vector of all data (population-level)
          n_obs=length(Y) # number total observations
     )
     
     # Population-level bayesian neg bin
     time_to_conf_model <- "
          model {
               
               for(i in 1:n_obs) {
                    
                    Y[i] ~ dnegbin(p, r)
               }
               
               # Priors
               p ~ dbeta(1, 1)
               r ~ dgamma(0.01, 0.01)
               
               mu <- (r*(1-p))/p
               sigma <- (r*(1-p))/(p^2)
          }
     "
     params <- c('p', 'r', 'mu', 'sigma')
     
     init.list <- replicate(n_chains, 
                            list(.RNG.name='lecuyer::RngStream',
                                 .RNG.seed= sample(1:1e6, 1)), 
                            simplify=FALSE)
     
     run.jags(model=time_to_conf_model,
              data=jags_data,
              monitor=params,
              n.chains=n_chains,
              adapt=1000,
              burnin=n_burn,
              sample=n_samp,
              thin=n_thin,
              inits=init.list,
              modules=c('lecuyer'),
              method="parallel",
              summarise=FALSE)
}



##' Fit spatial time to confirmation model
##' 
##' The function fits a hierarchical bayesian model to estimate the delay from time of symptom onset to time of case confirmation. The model
##' uses a negative binomial error structure with probability p and dispersion parameter r. Province-level distributions are also
##' fit with population-level hyper-paramters. This allows each province to diverge from population mean based on the data. Provinces with
##' little data regress to pop mean. Model fit to data using JAGS and 'runjags' pacakge. Defualt MCMC params should be increased for full run.
##' 
##' @param data data frame of linelist data produced by the 'read_MK_direct' function
##' @param n_chains number of chains to run (default=4)
##' @param n_burn number burn in samples (default=1000)
##' @param n_samps number of samples (default=1000)
##' @param n_thin interval to thin samples (default=1)
##' 
##' @return mcmc model object
##' 

fit_spatial_time_to_conf <- function(data,
                                     n_chains=4, # number chains
                                     n_burn=1000, # burn in
                                     n_samp=1000, # samples
                                     n_thin=1 # thin
) {
     
     require(runjags)
     
     tmp <- lapply(split(data$time_to_conf, data$province), function(x) x[!is.na(x)])
     Y_province <- do.call(rbind, lapply(tmp, "length<-", max(lengths(tmp))))
     Y <- as.vector(Y_province[!is.na(Y_province)])
     
     jags_data <- list(
          Y=Y, # vector of all data (population-level)
          Y_province=Y_province, # matrix of province-level data
          n_obs=length(Y), # number total observations 
          n_province=dim(Y_province)[1], # number of provinces
          n_obs_province=dim(Y_province)[2] # max observations in any province
     )
     
     time_to_conf_model <- "
          model {
               
               # Population-level model
               for (i in 1:n_obs) {
                    
                    Y[i] ~ dnegbin(p, r)
               }
               
               # Population-level hyper-priors
               p ~ dbeta(1, 1)
               r ~ dgamma(0.01, 0.01)
               
               mu <- (r*(1-p))/p # Population-level mean for time to confirmation
               
               # Province-level model
               for (j in 1:n_province) {
                    for (k in 1:n_obs_province) {
                         
                         Y_province[j,k] ~ dnegbin(p_province[j], r_province[j])
                    }
                    
                    # Province-level priors (with population-level hyper-parameters)
                    # Province-level estimates regress to population mean when province-level data are lacking
                    p_province[j] ~ dnorm(p, 100) T(0,1) 
                    r_province[j] ~ dnorm(r, 100) T(0,)
                    
                    mu_province[j] <- (r_province[j]*(1-p_province[j]))/p_province[j] # Province-level mean for time to confirmation
               }
          }
     "
     
     params <- c('p', 'r', 'mu', 'p_province', 'r_province', 'mu_province')
     
     init.list <- replicate(n_chains, 
                            list(.RNG.name='lecuyer::RngStream',
                                 .RNG.seed= sample(1:1e6, 1)), 
                            simplify=FALSE)
     
     run.jags(model=time_to_conf_model,
              data=jags_data,
              monitor=params,
              n.chains=n_chains,
              adapt=1000,
              burnin=n_burn,
              sample=n_samp,
              thin=n_thin,
              inits=init.list,
              modules=c('lecuyer'),
              method="parallel",
              summarise=FALSE)
}



##' Simulate spatial time to confirmation model
##' 
##' Function to simulate randome replicates of the time to confirmation estimated by the 'fit_spatial_tim_to_conf' function.
##' 
##' @param n number of replicates
##' @param mod model summary of posterior parameter estimates
##' @param province numeric indicator of which province to simulate values for (default=NULL returns all)
##' 
##' @return data frame where rows are province and columns are the n simulated values of time to confirmation
##' 

sim_spatial_time_to_conf <- function(n,
                                     mod,
                                     province=NULL
) {
     
     require(foreach)
     
     # Extract names of parameters and provinces
     tmp <- do.call(rbind, strsplit(row.names(mod), '\\[|\\]'))
     param_names <- tmp[,1]
     param_province <- suppressWarnings(as.numeric(tmp[,2]))
     
     # Reduce to focal province
     if (!is.null(province)) {
          
          mod <- mod[param_province %in% province,]
          tmp <- do.call(rbind, strsplit(row.names(mod), '\\[|\\]'))
          param_names <- tmp[,1]
          param_province <- suppressWarnings(as.numeric(tmp[,2]))
     }
     
     # n replicates foreach province
     out <- foreach(i=1:n, .combine='cbind') %do% {
          
          # Allow negbin params to vary within 95% confidence intervals
          pars <- cbind(
               truncnorm::rtruncnorm(n=1, 
                                     a=mod[param_names == 'p_province', 'Lower95'], 
                                     b=mod[param_names == 'p_province', 'Upper95'], 
                                     mean=mod[param_names == 'p_province', 'Mean'], 
                                     sd=mod[param_names == 'p_province', 'SD']),
               truncnorm::rtruncnorm(n=1, 
                                     a=mod[param_names == 'r_province', 'Lower95'], 
                                     b=mod[param_names == 'r_province', 'Upper95'], 
                                     mean=mod[param_names == 'r_province', 'Mean'], 
                                     sd=mod[param_names == 'r_province', 'SD'])
          )
          
          apply(pars, 1, function(x) rnbinom(1, prob=x[1], size=x[2]))
     }
     
     colnames(out) <- NULL
     data.frame(province=param_province[param_names == 'p_province'], out)
}


##' Simulate infection times
##' 
##' This function takes incidence data for one province and simulates n scenarios of augmented infection times for each case. 
##' Augmented infection times are gathered in a table and merged by date. Designed to by applied to a split data frame of incidence data.
##' 
##' @param n number of replicates
##' @param inc_data incidence data for one province
##' @param time_to_conf_model summary of posterior parameter estimates from spatial time to confirmation model
##' @param province_names vector of province names
##' @param n_cores number of cores to use when running n simulations in parallel (defualt=NULL which runs sequentially)
##' 
##' @return data frame with the province name, date, and additional columns containing each of the simulated scenarios
##' 
##' 
sim_infection_times <- function(n,
                                inc_data, 
                                time_to_conf_model,
                                province_names,
                                n_cores=NULL
) {
        
        if (length(unique(inc_data$Province_State)) > 1) stop('This function simulates infections times for one province at a time')
        
        if (!getDoParRegistered()) {
                
                if (!is.null(n_cores)) {
                        
                        clust <- parallel::makeCluster(n_cores)
                        
                        doParallel::registerDoParallel(clust)
                        parallel::clusterExport(clust, 
                                                c("sim_inf_period", "sim_spatial_time_to_conf", 
                                                  ls(environment())), envir=environment())
                }
                
                if (getDoParRegistered()) message(paste(getDoParWorkers(), 'cores registered for parallel backend')) 
        }
        
        message(paste0('Simulating ', n, ' scenarios of infection time in ', inc_data$Province_State[1]), '...')
        
        out <- foreach(j=1:n, .combine=merge_infection_time_sims) %dopar% { 
                
                as.Date(na.omit(do.call(
                        c, 
                        apply(inc_data, 1, function(inc_obs){
                                
                                n_cases <- round(as.numeric(inc_obs[names(inc_obs) == 'Incidence']))
                                
                                tryCatch({
                                        
                                        # Simulate incubation period for each case
                                        x <- sim_inf_period(n=n_cases)
                                        
                                        # Simulate time from symptom onset to case confirmation for each case
                                        y <- sim_spatial_time_to_conf(n=n_cases,
                                                                      mod=time_to_conf_model,
                                                                      province=which(province_names == inc_obs[names(inc_obs) == 'Province_State']))
                                        
                                        # Augmented infection times
                                        return(as.Date(inc_obs[names(inc_obs) == 'Date']) - (x + t(y[,-1])))
                                        
                                }, error = function (e) {
                                        
                                        return(as.Date(NA))
                                })
                        })
                )))
        }
        
        data.frame(province=as.character(inc_data$Province_State[1]),
                   date=dimnames(out)[[1]],
                   out,
                   row.names=NULL)
}

##' Merge infection time simulations
##' 
##' Helper function to merge raw simulated infection times into a matrix of infection incidence per day
##' 
##' @param sim1 first or previous simulations
##' @param sim2 new simulation to add
##' 
##' @return matrix with each column representing one simulation
##' 
##' 
merge_infection_time_sims <- function(sim1, sim2) {
     
     require(lubridate)
     
     if (is.Date(sim1)) sim1 <- table(sim1)
     if (is.Date(sim2)) sim2 <- table(sim2)
     
     tmp <- merge(as.matrix(sim1), as.matrix(sim2), by='row.names', all=T)
     out <- as.matrix(tmp[,-1])
     dimnames(out) <- list(tmp[,1])
     out
}





##' Simulate augmented infection incidence
##' 
##' This function simulates augmented infection times for each observed case in each province. 
##' 
##' @param n number of replicates
##' @param incidence_data incidence data for all provinces
##' @param linelist_data linelist data used to estimate time to confirmation
##' @param time_to_conf_model summary of posterior parameter estimates from spatial time to confirmation model (is NULL, fits this model)
##' @param summarise return the mean and 95% confidence intervals for infection incidence on each date (default = TRUE). If FALSE, returns all simulations
##' @param n_cores number of cores to use when running n simulations in parallel (defualt=NULL which runs sequentially)
##' 
##' @return data frame with the province name, date, and additional columns containing either mean and 95% CI or columns giving all simulations
##' 
##' 
augment_incidence <- function(n,
                              incidence_data,
                              linelist_data,
                              time_to_conf_model=NULL,
                              summarise=TRUE,
                              n_cores=NULL
) {
        
        if (is.null(time_to_conf_model)) {
                
                message('Fitting time to confirmation model...')
                time_to_conf_model <- summary(fit_spatial_time_to_conf(data=linelist_data,
                                                                       n_chains=4,
                                                                       n_burn=5000,
                                                                       n_samp=1000,
                                                                       n_thin=5))
        }
        
        if (!is.null(n_cores)) {
                
                clust <- parallel::makeCluster(n_cores)
                
                doParallel::registerDoParallel(clust)
                parallel::clusterExport(clust, 
                                        c("sim_inf_period", "sim_spatial_time_to_conf", 
                                          ls(environment())), envir=environment())
        }
        
        if (getDoParRegistered()) message(paste(getDoParWorkers(), 'cores registered for parallel backend')) 
        
        
        # Get provinces used in time_to_conf_model
        tmp <- lapply(split(linelist_data$time_to_conf, linelist_data$province), function(x) x[!is.na(x)])
        Y_province <- do.call(rbind, lapply(tmp, "length<-", max(lengths(tmp))))
        
        # Simulate infection times
        out <- do.call(rbind,
                       lapply(split(incidence_data_china, incidence_data_china$Province_State), 
                              function(x) sim_infection_times(n=n, 
                                                              inc_data=x, 
                                                              time_to_conf_model=time_to_conf_model,
                                                              province_names=dimnames(Y_province)[[1]],
                                                              n_cores=n_cores))
        )
        
        row.names(out) <- NULL
        colnames(out)[1:2] <- c('Province_State', 'Date')
        
        
        # Summarise simulations
        if (summarise == TRUE) {
                
                ci <- apply(out[,-c(1:2)], 1, function(x) quantile(x, probs=c(0.0275,0.975), na.rm=T))
                
                return(data.frame(out[,1:2],
                                  mean=apply(out[,-c(1:2)], 1, function(x) mean(x, na.rm=T)),
                                  lo95=ci[1,],
                                  hi95=ci[2,]))
        } else {
                
                return(out)
        }
}



##' Function to merge augmented infection incidence with original cumulative case counts
##' 
##' @param inf_data augmented infections times produced by augment_incidence function
##' @param cum_data a data frame with cumulative case data produced by the 'read_JHUCSSE_cases' function
##' @param inc_data a data frame produced by the 'est_daily_incidence' fucntion
##'
##' @return a data frame with both cumulative counts and inferred incidence
##'
merge_cum_inc <- function(inf_data,
                          cum_data,
                          inc_data
) {
     
     require(foreach)
     
     foreach(i=1:nrow(inc_data), .combine='rbind') %do% {
          
          sel <- which(format(cum_data$Update, '%Y-%m-%d') == format(inc_data$Date[i], '%Y-%m-%d') & cum_data$Province_State == inc_data$Province_State[i])
          tmp <- cum_data[sel,]
          
          if (nrow(tmp) == 0) {
               
               # Fill NAs for province/dates with no confirmed cases reported
               tmp <- cum_data[which(cum_data$Province_State == inc_data$Province_State[i])[1],]
               tmp[,3:7] <- NA
               
          } else if (nrow(tmp) > 1) {
               
               # If there are multiple entrys, use the confirmed counts with the latest time of day
               tmp <- tmp[order(tmp$Update, decreasing=TRUE)[1],]
          } 
          
          data.frame(tmp[,1:2],
                     Date=format(inc_data$Date[i], '%Y-%m-%d'),
                     tmp[,4:7],
                     Incidence=inc_data$Incidence[i],
                     row.names=NULL)
     }
}
