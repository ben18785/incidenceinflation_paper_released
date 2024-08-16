functions {
  real gamma_mean_sd_lpdf(real x, real mu, real sigma) {
    return gamma_lpdf(x | mu ^ 2 / (sigma ^ 2), mu / (sigma ^ 2));
  }

  real gamma_cdf_mean_variance(real x, real mu, real sigma) {
    return gamma_cdf(x, mu ^ 2 / (sigma ^ 2), mu / (sigma ^ 2));
  }

  real p_detect(real t, real t_prime, real t_double_prime, row_vector theta,
                int is_gamma_reporting_delay) {
    real S_t_prime;
    real S_t_double_prime;
    real loc = theta[1];
    real scale = theta[2];

    if (is_gamma_reporting_delay == 1) {
      S_t_prime = 1 - gamma_cdf_mean_variance(t_prime - t, loc, scale) + 1e-10; // add a bit to avoid roundoff
      S_t_double_prime = 1 - gamma_cdf_mean_variance(t_double_prime - t, loc, scale);
    } else if (is_gamma_reporting_delay == 0) {
      S_t_prime = 1 - lognormal_cdf(t_prime - t, loc, scale);
      S_t_double_prime = 1 - lognormal_cdf(t_double_prime - t, loc, scale);
    }

    return (S_t_prime - S_t_double_prime) / S_t_prime;
  }
}

data {
  int N_delay;
  array[N_delay] real time_reported;
  array[N_delay] real time_onset;
  array[N_delay] int cases_reported;
  array[N_delay] int cases_true;
  int n_reporting_window;
  array[N_delay] int<lower=1, upper=n_reporting_window> reporting_window;
  
  //options
  int<lower=0, upper=1> is_gamma_reporting_delay; // 1-> gamma; 0-> lognormal
  int<lower=0, upper=1> is_binomial; // 1-> binomial; 0-> beta-binomial
  
  int N_sim;
  vector[N_sim] reporting_days_sim;
}

parameters {
  matrix<lower=0>[n_reporting_window, 2] theta;
  array[is_binomial==1 ? 0 : 1] real<lower=1> kappa;
}

model {
  // likelihood
  for (t in 1 : (N_delay - 1)) {
      if (time_onset[t] == time_onset[t + 1]) {
        int cases_observed = cases_reported[t + 1] - cases_reported[t];
        int cases_remaining = cases_true[t] - cases_reported[t];
        row_vector[2] theta_t = theta[reporting_window[t]];
        real p = p_detect(time_onset[t], time_reported[t] + 0.5,
                          time_reported[t + 1] + 0.5, theta_t,
                          is_gamma_reporting_delay);
        if(is_binomial)
          cases_observed ~ binomial(cases_remaining, p);
        else
          cases_observed ~ beta_binomial(cases_remaining, p * kappa[1], (1 - p) * kappa[1]);
      }
  }
  if(is_binomial == 0)
    kappa ~ cauchy(0, 10);
  for(i in 1:n_reporting_window)
    theta[i] ~ cauchy(5, 10);
}

generated quantities {
  matrix[N_sim, n_reporting_window] m_ecdf_delay;
  vector[N_delay - 1] v_ecdf;
  vector[N_delay - 1] v_cases_observed;
  vector[N_delay - 1] v_cases_remaining;
  vector[N_delay] v_cases_reported;
  vector[N_delay - 1] v_cases_reported_sim;
  vector[N_delay - 1] v_p;
  for(i in 1:n_reporting_window) {
    for(j in 1:N_sim) {
      row_vector[2] theta_t = theta[i];
      if(is_gamma_reporting_delay == 1)
        m_ecdf_delay[j, i] = gamma_cdf_mean_variance(reporting_days_sim[j], theta_t[1], theta_t[2]);
      else
        m_ecdf_delay[j, i] = lognormal_cdf(reporting_days_sim[j], theta_t[1], theta_t[2]);
    }
  }
  
  for (t in 1 : (N_delay - 1)) {
    if (time_onset[t] == time_onset[t + 1]) {
        int cases_observed = cases_reported[t + 1] - cases_reported[t];
        int cases_remaining = cases_true[t] - cases_reported[t];
        row_vector[2] theta_t = theta[reporting_window[t]];
        real p = p_detect(time_onset[t], time_reported[t] + 0.5,
                          time_reported[t + 1] + 0.5, theta_t,
                          is_gamma_reporting_delay);
        v_cases_observed[t] = cases_observed;
        v_cases_remaining[t] = cases_remaining;
        v_p[t] = p;
   if(is_binomial)
        v_ecdf[t] = binomial_cdf(cases_observed|cases_remaining, p);
   else
        v_ecdf[t] = beta_binomial_cdf(cases_observed|cases_remaining, p * kappa[1], (1 - p) * kappa[1]);
    }
  }
  
  int cases_remaining_tmp, cases_observed_tmp;
  for(t in 1: (N_delay - 1)) {
    if (t == 1) {
      cases_remaining_tmp = cases_true[t];
      v_cases_reported[t] = 0;
    } else if (time_onset[t] != time_onset[t - 1]) {
      cases_remaining_tmp = cases_true[t];
      v_cases_reported[t] = 0;
    }
    
    if (time_onset[t] == time_onset[t + 1]) {
      row_vector[2] theta_t = theta[reporting_window[t]];
      real p = p_detect(time_onset[t], time_reported[t] + 0.5,
                          time_reported[t + 1] + 0.5, theta_t,
                          is_gamma_reporting_delay);
      if(is_binomial)
        cases_observed_tmp = binomial_rng(cases_remaining_tmp, p);
      else
        cases_observed_tmp = beta_binomial_rng(cases_remaining_tmp, p * kappa[1], (1 - p) * kappa[1]);
        
      cases_remaining_tmp -= cases_observed_tmp;
      v_cases_reported_sim[t] = cases_observed_tmp;
      v_cases_reported[t + 1] = v_cases_reported[t] + cases_observed_tmp;
    }
  }
}
