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
    int<lower=0> test_pos[N];
    vector<lower=0>[N] t_symp_test;
    int<lower=1> T_max;
    int<lower=1> exposed_n;
    int<lower=0> exposed_pos;
    // vector<lower=0,upper=T_max>[N] t_min;
    // vector<lower=0,upper=T_max>[N] t_max;
}

// 'z_survey_pos' is the total number of positive predictions for the survey.
// We also find the log-time since infection for the cohort 't_log', with
// its orthogonal squared and cubic terms ('t_log2' and 't_log3'). We will
// later estimate the sensitivity for all values from t=1:T_max, for which
// we need 't_log_survey', 'tls_2', and 'tls_3'.
transformed data {
    vector[T_max] t_new;
    real t_mean;
    real t_sd;

    for(i in 1:T_max){
        t_new[i] = log(i);
    }
    t_mean=sum(t_new)/T_max;
    t_sd=sd(t_new);
}

// 'p' is our parameter of interest, the true seroincidence. 'spec' is the
// specificity of the RF predictions. The 'beta' values are the coefficients
// of the logistic regression model for finding the time-varying sensitivity
// (below). 'd' is the probability of infection 1:T_max days ago.
parameters {
    real beta_0;
    real beta_1;
    real beta_2;
    real beta_3;
    real<lower=0, upper=1> attack_rate;
    vector<lower=0>[N] t_exp_symp;
}

// 'sens_tv' is the time-varying sensitivity at each 'd' from 1:T_max. 'sens' is
// the overall sensitivity integrated across the estimated days back 'd'.
transformed parameters{
    vector[N] t;
    vector[N] t_ort;
    vector[N] t_ort2;
    vector[N] t_ort3;

    t = t_symp_test + t_exp_symp;
    for(i in 1:N){
        t_ort[i] = (log(t[i])-t_mean)/t_sd;
        t_ort2[i] = t_ort[i]^2;
        t_ort3[i] = t_ort[i]^3;
    }

}

//  We observe 'z_survey' cases as a binomial distribution based on the
//  survey sample size with observations coming as the sum of the true
//  positive rate p*sens and false negative rate (1-p)*(1-spec). Sensitivity
//  is modeled as a logistic regression with a cubic polynomial for log-time.
//  We assume that the specificity is distributed as a binomial.
model {
    exposed_pos ~ binomial(exposed_n, attack_rate);
    test_pos ~ bernoulli_logit(beta_0+beta_1*t_ort+beta_2*t_ort2+beta_3*t_ort3);

    t_exp_symp ~ gamma(5.81, 1/0.95);
}

generated quantities{
    vector<lower=0, upper=1>[T_max] sens;
    vector<lower=0, upper=1>[T_max] npv;
    vector[T_max] t_new_ort;
    vector[T_max] t_new_ort2;
    vector[T_max] t_new_ort3;
    vector[N] log_lik;

    t_new_ort = (t_new-t_mean)/t_sd;
    for(i in 1:T_max){
        t_new_ort2[i] = t_new_ort[i]^2;
        t_new_ort3[i] = t_new_ort[i]^3;
    }

    // sens=inv_logit(beta_0+beta_1*t_new_ort+beta_2*t_new_ort2);
    sens=inv_logit(beta_0+beta_1*t_new_ort+beta_2*t_new_ort2+beta_3*t_new_ort3);
    for(i in 1:T_max){
        npv[i]=(1-attack_rate)/((1-sens[i])*attack_rate+(1-attack_rate));
    }

    for(i in 1:N){
        log_lik[i] = bernoulli_logit_lpmf(test_pos[i] | beta_0+beta_1*t_ort[i]+beta_2*t_ort2[i]+beta_3*t_ort3[i]);
    }
}
