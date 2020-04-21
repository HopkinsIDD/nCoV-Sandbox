//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic.
// We assume that the sensitivity of the diagnostic decays with time.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// We have input data from a cohort (which determines the diagnostic) and
// the survey (which we'd like to estimate). 'z_survey' are the predictions
// of the survey seroincidence from a random forest model that has 'N_survey'
// observations. 'N_cohort_pos' are the number of seroincident observations
// in the cohort data. 't' is the time since infection for each seropositive
// cohort observation and 'z_cohort_pos' are the predictions of these
// observations. There were 'N_cohort_neg' observations of seronegatives in the
// cohort study, which were predicted as 'z_cohort_neg'.
data {
    int<lower=1> N;
    int<lower=1> T_max;
    int<lower=0> test_pos[N];
    int<lower=1> test_n[N];
    int<lower=0> t_symp_test[N];
    int<lower=1> exposed_n;
    int<lower=0> exposed_pos;
}

// 'z_survey_pos' is the total number of positive predictions for the survey.
// We also find the log-time since infection for the cohort 't_log', with
// its orthogonal squared and cubic terms ('t_log2' and 't_log3'). We will
// later estimate the sensitivity for all values from t=1:T_max, for which
// we need 't_log_survey', 'tls_2', and 'tls_3'.
transformed data {
    real t[N];
    real t_mean;
    real t_sd;
    real t_ort[N];
    real t_new[T_max];
    real t_new_ort[T_max];

    for(i in 1:T_max)
        t_new[i] = log(i);

    t_mean = sum(t_new)/T_max;
    t_sd = sd(t_new);

    for(i in 1:N){
        t[i] = t_symp_test[i]+5;
        t_ort[i] = (log(t[i])-t_mean)/t_sd;
    }
    for(i in 1:(T_max)){
        t_new_ort[i] = (log(i)-t_mean)/t_sd;
    }

}

// 'p' is our parameter of interest, the true seroincidence. 'spec' is the
// specificity of the RF predictions. The 'beta' values are the coefficients
// of the logistic regression model for finding the time-varying sensitivity
// (below). 'd' is the probability of infection 1:T_max days ago.
parameters{
    real beta_0;
    real beta_1;
    real beta_2;
    real beta_3;
    real<lower=0, upper=1> attack_rate;
}

transformed parameters{
    real<lower=0> db_dt[3];

    for(i in 1:3)
        db_dt[i] = beta_1+2*beta_2*(log(i+1)-t_mean)/t_sd+3*beta_3*((log(i+1)-t_mean)/t_sd)^2;
}

//  We observe 'z_survey' cases as a binomial distribution based on the
//  survey sample size with observations coming as the sum of the true
//  positive rate p*sens and false negative rate (1-p)*(1-spec). Sensitivity
//  is modeled as a logistic regression with a cubic polynomial for log-time.
//  We assume that the specificity is distributed as a binomial.
model {
    exposed_pos ~ binomial(exposed_n, attack_rate);
    for(i in 1:N){
        test_pos[i] ~ binomial_logit(test_n[i], beta_0+beta_1*t_ort[i]+beta_2*t_ort[i]^2+beta_3*t_ort[i]^3);
    }
    // beta_0 ~ normal(0,1);
    // beta_1 ~ normal(0,1);
    // beta_2 ~ normal(0,1);
    // beta_3 ~ normal(0,1);
}

generated quantities{
    vector<lower=0, upper=1>[T_max] sens;
    vector<lower=0, upper=1>[T_max] npv;
    vector[N] log_lik;

    for(i in 1:T_max){
        sens[i] = inv_logit(beta_0+beta_1*t_new_ort[i]+beta_2*t_new_ort[i]^2+beta_3*t_new_ort[i]^3);
    }

    for(i in 1:T_max){
        npv[i]=(1-attack_rate)/((1-sens[i])*attack_rate+(1-attack_rate));
    }

    for(i in 1:N){
        log_lik[i] = binomial_logit_lpmf(test_pos[i] | test_n[i], beta_0+beta_1*t_ort[i]+beta_2*t_ort[i]^2+beta_3*t_ort[i]^3);
    }
}
