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
##' @param combine_CIs should we combine the confidence intervals
##'    into a single entry, or hav two separate columns
##'
##' @return a data frame with includein colums with OR and CI
##'
OR_table_age <- function(data, combine_CIs=TRUE) {
  age_odds <- data %>%
    drop_na(age_cat)%>%
    group_by(age_cat,death)%>%
    summarize(n())%>%
    mutate(death=ifelse(death,"dead","alive"))%>%
    pivot_wider(names_from = death,values_from=`n()`) %>%
    replace_na(list(dead=0))

  if ("(50,60]"%in%data$age_cat) {
    data$age_cat <- relevel(data$age_cat, ref="(50,60]")
  } else {
    data$age_cat <- relevel(data$age_cat, ref="50-59")
  }

  age_odds_mdl <- glm(death~age_cat, data=data,
                      family=binomial)

  n_coefs <- length(age_odds_mdl$coef)
  age_odds$OR <- c(exp(age_odds_mdl$coef[2:6]),1,
                   exp(age_odds_mdl$coef[7:n_coefs]))

  if (combine_CIs) {
    tmp <- round(exp(confint(age_odds_mdl)),2)
    tmp <- apply(tmp, 1, paste, collapse=",")
    age_odds$CI <- c(tmp[2:6],"-",tmp[7:n_coefs])
  } else {
      tmp <-exp(confint(age_odds_mdl))
      age_odds$CI_low = c(tmp[2:6,1],1,tmp[7:n_coefs,1])
      age_odds$CI_high = c(tmp[2:6,2],1,tmp[7:n_coefs,2])
  }

  return(age_odds)

}


##'
##' Function to extract approximate epidemic curves
##' from the cumulative case data.
##'
##' @param cum_data a data frame with cumulative case data in oit
##' @param first_date the first date to infer
##' @param second_date the last date to infer...shold be iwthin data range
##'
##'
##' @return a data frame with roughly estimated incidence in it
##'
est_daily_incidence <- function (cum_data,
                                 first_date,
                                 last_date) {
  analyze <-   cum_data %>% drop_na(Confirmed)

  ##Get the implied daily incidence for each province
  ##by fitting a monitonically increasing spline and then
  ##taking the difference (a little less sensitive to
  ##perturbations in reporting than taking raw difference).
  ##Making sure only to infer over trhe suport
  tmp_dt_seq <- seq(first_date, last_date, "days")
  incidence_data<- analyze %>% nest(-Province_State) %>%
    mutate(cs=map(data, ~splinefun(x=.$Update, y=.$Confirmed,
                                   method="hyman"))) %>%
    mutate(Incidence=map2(cs,data, ~data.frame(Date=tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)],
                                               Incidence= diff(c(0, pmax(0,.x(tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)]))))))) %>%
    unnest(Incidence) %>% select(-data) %>% select(-cs)

  return(incidence_data)

}


##'
##' Function to extract approximate epidemic curves
##' from the cumulative case data.
##'
##' @param cum_data a data frame with cumulative case data in oit
##' @param first_date the first date to infer
##' @param second_date the last date to infer...shold be iwthin data range
##'
##'
##' @return a data frame with roughly estimated incidence in it
##'
est_daily_deaths <- function (cum_data,
                                 first_date,
                                 last_date) {
  analyze <-   cum_data %>% drop_na(Deaths)

  ##Get the implied daily incidence for each province
  ##by fitting a monitonically increasing spline and then
  ##taking the difference (a little less sensitive to
  ##perturbations in reporting than taking raw difference).
  ##Making sure only to infer over trhe suport
  tmp_dt_seq <- seq(first_date, last_date, "days")
  death_data<- analyze %>% nest(-Province_State) %>%
    mutate(cs=map(data, ~splinefun(x=.$Update, y=.$Deaths,
                                   method="hyman"))) %>%
    mutate(Deaths=map2(cs,data, ~data.frame(Date=tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)],
                                               Deaths= diff(c(0, pmax(0,.x(tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)]))))))) %>%
    unnest(Deaths) %>% select(-data) %>% select(-cs)

  return(death_data)

}




##'
##' Function to extract approximate epidemic curves
##' from the cumulative case data.
##'
##' @param cum_data a data frame with cumulative case data in oit
##' @param first_date the first date to infer
##' @param last_date the last date to infer...shold be iwthin data range
##' @param na_to_zeros
##'
##'
##' @return a data frame with roughly estimated incidence in it
##'
est_daily_recovered <- function (cum_data,
                              first_date,
                              last_date,
                              na_to_zeros=FALSE) {
    if (na_to_zeros) {
        analyze <-   cum_data %>% replace(is.na(.), 0)
    } else {
        analyze <-   cum_data %>% drop_na(Recovered)
    }

  ##Get the implied daily incidence for each province
  ##by fitting a monitonically increasing spline and then
  ##taking the difference (a little less sensitive to
  ##perturbations in reporting than taking raw difference).
  ##Making sure only to infer over trhe suport
  tmp_dt_seq <- seq(first_date, last_date, "days")
  recovered_data<- analyze %>% nest(-Province_State) %>%
    mutate(cs=map(data, ~splinefun(x=.$Update, y=.$Recovered,
                                   method="hyman"))) %>%
    mutate(Incidence=map2(cs,data, ~data.frame(Date=tmp_dt_seq[tmp_dt_seq>=min(.y$Update)],
                                               Recovered= diff(c(0, pmax(0,.x(tmp_dt_seq[tmp_dt_seq>=min(.y$Update)]))))))) %>%
    unnest(Incidence) %>% select(-data) %>% select(-cs)

  return(recovered_data)

}

