##' Make a ggplot object that shows the age distirbution of cases
##' with indicator for living or dead.
##' 
##' @param data the nCoV data to run on. Assued to have age categories.
##' 
##' @return a ggplot object showing the age distribution of cases
##' 
age_dist_graph <- function (data) {
  require(ggplot2)
  rc <- ggplot(drop_na(data, age_cat), 
           aes(x=age_cat, fill=as.factor(death))) + 
    geom_bar( color="grey") + coord_flip() + xlab("Age Catergory")
  
  return(rc)
  
}

##' Make the epi curve of reported cases. With living and 
##' dead indicated.
##' 
##' @param data the nCov data to run on
##' 
##' @return a ggplot object showing the epidemiologic curve.
##'
epi_curve_reported_cases <- function(data) {
  ggplot(data, aes(x=symptom_onset, fill=as.factor(death))) +
    geom_bar()
}


##'
##' Make a table capturing the odds ratio of death by geneder
##' 
##' @param data
##' 
##' @return a data frame with includein colums with OR and CI
##' 
OR_table_gender <- function(data) {
  gender_odds <- data%>%group_by(gender, death)%>%
    summarize(n())%>%
    mutate(death=ifelse(death,"dead","alive"))%>%
    pivot_wider(names_from = death,values_from=`n()`)
  
  #%>%
  #mutate(OR=dead/alive)
  
  #gender_odds <- gender_odds%>%mutate(OR=OR/gender_odds$OR[1])
  gender_odds_mdl <- glm(death~gender, data=data,
                         family=binomial)
  gender_odds$OR <- c(1, exp(gender_odds_mdl$coef[2]))
  gender_odds$CI <- c("-",
                      paste(round(exp(confint(gender_odds_mdl)[,2]),2),
                            collapse = ","))
  return(gender_odds)
}


##'
##' Make a table capturing the odds ratio of death by age cat
##' 
##' @param data the data inclding an age_cat column
##' 
##' @return a data frame with includein colums with OR and CI
##' 
OR_table_age <- function(data) {
  age_odds <- data %>%
    drop_na(age_cat)%>%
    group_by(age_cat,death)%>%
    summarize(n())%>%
    mutate(death=ifelse(death,"dead","alive"))%>%
    pivot_wider(names_from = death,values_from=`n()`) %>%
    replace_na(list(dead=0))
  
  data$age_cat <- relevel(data$age_cat, ref="(50,60]")
  age_odds_mdl <- glm(death~age_cat, data=data,
                      family=binomial)
  
  age_odds$OR <- c(exp(age_odds_mdl$coef[2:6]),1,
                   exp(age_odds_mdl$coef[7:9]))

  tmp <- round(exp(confint(age_odds_mdl)),2)
  tmp <- apply(tmp, 1, paste, collapse=",")
  age_odds$CI <- c(tmp[2:6],"-",tmp[7:9])
  
  return(age_odds)
  
}