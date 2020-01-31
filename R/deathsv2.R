### get survival curve
library(survival)
library(ggplot2)
setwd("/Users/jkedwar/Box Sync/Research/epidemics/ncov19/ncov19git/nCoV-Sandbox")


#' A function to filter
lim <- function(df){
  df <- df[df$Country %in% c("China", "Mainland China"),]
  return(df)
}

#' A function to calculate new deaths in this time step
#' @param df1 dataframe with counts of deaths occuring up until today
#' @param df2 data frame with counts of deaths occuring up until tomorrow
#' @return dataframe containing death that occurred today
# assume constant hazard within some window (can put shape on later)
#this function outputs probability of dying given that you are still alive on this date
getdeltadeaths <- function(df1, df2){
  datint <- merge(df1, df2, by = "Province", all = TRUE)
  #datint$deltacases <- datint$Confirmed.y - ifelse(is.na(datint$Confirmed.x), 0, datint$Confirmed.x)
  datint$deltadeaths <- datint$Deaths.y - ifelse(is.na(datint$Deaths.x), 0, datint$Deaths.x)
  datint$alivecases <- datint$Confirmed.x - ifelse(is.na(datint$Deaths.x), 0, datint$Deaths.x)
  datint$pdeath <- datint$deltadeaths/datint$alivecases
  datint$pdeath <- ifelse(is.na(datint$pdeath), 0, datint$pdeath)
  datint$Date <- as.Date(substr(as.character(datint$Date.x), 1, 9), format="%m/%d/%y")
  datint <- datint[!is.na(datint$Date.y), c(1, 13, 12)]
}


#' A function to calculate new cases in this time step
#' @param df1 dataframe with counts of deaths occuring up until yesterday
#' @param df2 data frame with counts of deaths occuring up until today
#' @return dataframe containing counts of new cases today
getdeltacases <- function(df1, df2){
  datint <- merge(dat1, dat2, by = "Province", all = TRUE)
  datint$deltacases <- datint$Confirmed.y - ifelse(is.na(datint$Confirmed.x), 0, datint$Confirmed.x)
  datint$deltacases <- ifelse(is.na(datint$deltacases), 0, datint$deltacases)
  #datint$deltadeaths <- datint$Deaths.y - ifelse(is.na(datint$Deaths.x), 0, datint$Deaths.x)
  datint <- datint[!is.na(datint$Date.y), c(1, 7, 10)]
}

#' A function to make a pseudo line list out of tabular data
#' @param df dataframe with counts of new cases occuring today
#' @return dataframe with one record per new case occuring today
getll <- function(df, pconfirmed = 1){
  datint <- df[rep(seq_len(nrow(df)), df$deltacases * 1/pconfirmed),  ]
  datint$date <- as.Date(substr(as.character(datint[,2]), 1, 9), format="%m/%d/%y")
  datint<- datint[, c(1, 4)]
  datint$deathDate <- NaN
  datint$death <- NaN
  datint$confirmed <- rbinom(nrow(datint), 1, pconfirmed)
  names(datint) <- c("Province",  "onsetDate", "deathDate", "death", "confirmed")
  return(datint)
  
}

#' A function to impute random dates

rdate <- function(x,
                  min = paste0(format(Sys.Date(), '%Y'), '-01-01'),
                  max = paste0(format(Sys.Date(), '%Y'), '-12-31'),
                  sort = TRUE) {
  
  dates <- sample(seq(as.Date(min), as.Date(max), by = "day"), x, replace = TRUE)
  if (sort == TRUE) {
    sort(dates)
  } else {
    dates
  }
  
}


  
#' A function to simulate deaths occuring in line list
#' @param olddeaths dataframe with pseudo line list with dates of death up until yesterday
#' @param newdeaths dataframe with counts of new deaths today by province
#' @return dataframe with one record per case for the entire epidemic to date, with simulated dates of death, if applicable
adddeaths <- function(olddeaths, newdeaths){
  olddeaths <- olddeaths[, !names(olddeaths) %in% c("pdeath", "Date")]
  olddeaths$death <-  ifelse(is.na(olddeaths$death), 0, olddeaths$death)
  newdeathsll <- merge(olddeaths, newdeaths, by = c("Province"), all = TRUE)
  newdeaths$pdeath <-  ifelse(is.na(newdeaths$pdeath), 0, newdeaths$pdeath)
  newdeathsll$death <- ifelse(newdeathsll$death == 0 & newdeathsll$onsetDate <= newdeathsll$Date & newdeathsll$confirmed == 1, #only confirmed cases can die
                              rbinom(nrow(newdeathsll), 1, newdeathsll$pdeath), newdeathsll$death) # we could modify this line to allow hazard to vary over time since symptom onset
  # regarding the above line, need to think about whether we could use some sort of iterative algorithm to guess the
  # appropriate shape of the hazard function...
  newdeathsll$deathDate <- ifelse(newdeathsll$death == 1, newdeathsll$Date, NA)
  return(newdeathsll[!is.na(newdeathsll$onsetDate), c("Province", "onsetDate", "deathDate", "death", "confirmed")])
}

#' A function to update pseduo line list with deaths with new data each day
#' @param current dataframe with pseudo line list with dates of death up until yesterday
#' @param bfdat raw counts dataframe from yesterday (cumulative confirmed cases and cumulative deaths)
#' @param nowdat raw counts dataframe from today (cumulative confirmed cases and cumulative deaths)
#' @param nextdat raw counts dataframe from tomorrow (cumulative confirmed cases and cumulative deaths)
#' @return dataframe with one record per case for the entire epidemic to date, with simulated dates of death, if applicable
updatedat <- function(current, bfdat, nowdat, nextdat, pconfirmed = 1){
  datb <- getdeltacases(bfdat, nowdat)
  datc <- getdeltadeaths(nowdat, nextdat)
  datll <- getll(datb, pconfirmed = pconfirmed)
  linelist <- rbind(current, datll)
  newdeaths <- adddeaths(linelist, datc)
}

# process initiatial data as of 1/30
dat1 <- read.csv(file = "data/JHUCSSE Total Cases 1-22-2020.csv")[, c(1:4)]
dat2 <- read.csv(file = "data/JHUCSSE Total Cases 1-23-2020.csv")[, c(1:4, 7)]
dat3 <- read.csv(file = "data/JHUCSSE Total Cases 1-24-2020.csv")[, c(1:4, 7)]
dat4 <- read.csv(file = "data/JHUCSSE Total Cases 1-25-2020.csv")[, c(1:4, 7)]
dat5 <- read.csv(file = "data/JHUCSSE Total Cases 1-26-2020.csv")[, c(1:5)]
dat6 <- read.csv(file = "data/JHUCSSE Total Cases 1-27-2020.csv")[, c(1:5)]
dat7 <- read.csv(file = "data/JHUCSSE Total Cases 1-28-2020.csv")[, c(1:5)]

dat1$deaths <- 0


dats <- c("dat1", "dat2", "dat3", "dat4", "dat5", "dat6", "dat7")
for(df in dats){
  assign(df, setNames(get(df),  c("Province", "Country", "Date", "Confirmed", "Deaths")))
}

#limit to mainland china
dat1 <- lim(dat1)
dat2 <- lim(dat2)
dat3 <- lim(dat3)
dat4 <- lim(dat4)
dat5 <- lim(dat5)
dat6 <- lim(dat6)
dat7 <- lim(dat7)

#input confirmation rates (vector of length 6, here) - assuming 1s for now?
pconfirmed <- c(1, 1, 1, 1, 1, 1)
#pconfirmed <- c(.5, .25, .2, .2, .1, .1)

#start loop here
res <- matrix(nrow = 200, ncol = 1)
for(i in 1:200){
  #organize first file, which has different vars and slightly different structure than others
  dat1b <- dat1
  dat1b$Date.y <- dat1b$Date
  dat1b$deltacases <- dat1$Confirmed
  dat1b$deltacases <- ifelse(is.na(dat1b$deltacases), 0, dat1b$deltacases)
  dat1b <- dat1b[, c(1, 3, 7)]
  dat1c <- getdeltadeaths(dat1, dat2)
  dat1ll <- getll(dat1b, pconfirmed = pconfirmed[1])
  #assign random onset date (between 1/1/20 and 1/21/20) to prevalent cases on 1/21 (could alter this to make fewer cases early and more later...)
  dat1ll$onsetDate <- rdate(length(dat1ll$onsetDate), min = '2020-01-01', max = "2020-01-21", sort = FALSE)
  deaths1 <- merge(dat1ll, dat1c, by = c("Province"), all = TRUE)
  deaths1$death <- ifelse(deaths1$onsetDate<=deaths1$Date, rbinom(nrow(deaths1), 1, deaths1$pdeath), 0)
  deaths1$deathDate <- ifelse(deaths1$death == 1, deaths1$Date, "NA")
  deaths1 <- deaths1[!is.na(deaths1$onsetDate), c("Province", "onsetDate", "deathDate", "death", "confirmed")]
  
  # use function to update with new data each day
  deaths2 <- updatedat(deaths1, dat1, dat2, dat3, pconfirmed = pconfirmed[2])
  deaths3 <- updatedat(deaths2, dat2, dat3, dat4, pconfirmed = pconfirmed[3])
  deaths4 <- updatedat(deaths3, dat3, dat4, dat5, pconfirmed = pconfirmed[4])
  deaths5 <- updatedat(deaths4, dat4, dat5, dat6, pconfirmed = pconfirmed[5])
  deaths6 <- updatedat(deaths5, dat5, dat6, dat7, pconfirmed = pconfirmed[6])
  
  uptodate <- deaths6 #update when new data become avail
  
  # turn this on as a check on the number of deaths - right now we seemt to be returning the observed numbers of deaths
  # sum(deaths1$death, na.rm = TRUE)
  # sum(deaths2$death, na.rm = TRUE)
  # sum(deaths3$death, na.rm = TRUE)
  # sum(deaths4$death, na.rm = TRUE)
  # sum(deaths5$death, na.rm = TRUE)
  # sum(deaths6$death, na.rm = TRUE)
  
  uptodate$delta <- ifelse(is.na(uptodate$deathDate), 0, 1)
  uptodate$deathDate2 <- as.Date(as.numeric(uptodate$deathDate), origin = "1970-01-01")
  uptodate$t <- uptodate$deathDate2- uptodate$onsetDate #time to death
  uptodate$t <- ifelse(is.na(uptodate$t), Sys.Date() - uptodate$onsetDate, uptodate$t) #administratively censor today
  deathci <- survfit(Surv(t, delta == 1) ~ 1, data = uptodate)
  deathci2 <- data.frame( t = deathci$time, s = deathci$surv)
  deathci2$r <- 1-deathci2$s
  
  res[i, ] <- unlist(tail(deathci2[deathci2$t<=20,"r"], n=1)) #risk at 20 days after symptom onset (risk sets get small after, at the moment. will increase as time passes)
}

summary(res)
#how can we calibrate these results to account for the fact that we are using on confirmed (ie most severe) cases?

plot(hist(res, breaks = 50))

#plot the last draw for fun

deathpl <- ggplot()+
  geom_step(data = deathci2, aes(x = t, y = r), size = 1) +
  scale_y_continuous( name = "Proportion dying prior to time t") +
  coord_cartesian(xlim = c(0, 45), ylim = c(0, .5)) +
  scale_x_continuous( name = "Days after symptom onset") +
  theme(text = element_text(size = 15), legend.title = element_blank(), legend.text.align = 0, legend.position = c(.7, .85))
#,
deathpl

