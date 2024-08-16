functions {
  real gamma_mean_sd_lpdf(real x, real mu, real sigma) {
    return gamma_lpdf(x | mu ^ 2 / (sigma ^ 2), mu / (sigma ^ 2));
  }

  real gamma_cdf_mean_variance(real x, real mu, real sigma) {
    return gamma_cdf(x | mu ^ 2 / (sigma ^ 2), mu / (sigma ^ 2));
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
      S_t_prime = 1 - lognormal_cdf(t_prime - t | loc, scale);
      S_t_double_prime = 1 - lognormal_cdf(t_double_prime - t | loc, scale);
    }

    return (S_t_prime - S_t_double_prime) / S_t_prime;
  }
}
data {
  // renewal model data
  int N; // number of data points
  int K; // number of R segments
  array[N] int window; // assigns a time point to a given segment
  array[N] int C; // case series
  int wmax; // max day of generation time distribution
  row_vector[wmax] w; // generation distribution

  // delay data
  int N_delay;
  array[N_delay] real time_reported;
  array[N_delay] real time_onset;
  array[N_delay] int cases_reported;
  array[N_delay] int cases_true; // note this contains similar values to C except has repeats
  int n_reporting_window;
  array[N_delay] int<lower=1, upper=n_reporting_window> reporting_window;

  // options
  int<lower=0, upper=1> is_poisson; // 1-> Poisson; 0-> negative binomial
  int<lower=0, upper=1> is_gamma_reporting_delay; // 1-> gamma; 0-> lognormal
  int<lower=1, upper=2> prior_R_choice; // 1-> random walk; 2-> gamma
  real prior_R_a;
  real prior_R_b;
  real prior_kappa_a;
  real prior_kappa_b;
  real prior_theta_1_a;
  real prior_theta_1_b;
  real prior_theta_2_a;
  real prior_theta_2_b;
}
transformed data {
  int is_random_walk = prior_R_choice == 1 ? 1 : 0;
}
parameters {
  array[K] real<lower=0> R;
  matrix<lower=0>[n_reporting_window, 2] theta;
  array[is_poisson==1 ? 0 : 1] real<lower=0> kappa;
  array[is_random_walk==1 ? 1 : 0] real<lower=0> sigma;
}
transformed parameters {
  array[N - 1] real E_cases;

  {
    vector[wmax] I_temp;
    array[N - 1] real sum_w;
    for (t in 1 : N) {
      if (t > 1) {
        if (t <= wmax) {
          int kk = wmax - t + 1;
          for (i in 1 : (t - 1))
            I_temp[i] = C[t - i]; // needs to be lagged by one time point
          for (i in 1 : kk)
            I_temp[i + t - 1] = 0;
        } else {
          for (i in 1 : wmax)
            I_temp[i] = C[t - i]; // needs to be lagged by one time point
        }
        sum_w[t - 1] = w * I_temp;
        E_cases[t - 1] = R[window[t]] * sum_w[t - 1];
      }
    }
  }
}
model {
  // renewal model
  // likelihood
  if (is_poisson == 0) {
    C[2 : N] ~ neg_binomial_2(E_cases, rep_vector(kappa[1], N - 1));
  } else if (is_poisson == 1) {
    C[2 : N] ~ poisson(E_cases);
  }

  // priors
  if (is_poisson == 0)
    kappa ~ cauchy(prior_kappa_a, prior_kappa_b);

  if (prior_R_choice == 1) {
    // forward random walk

    sigma ~ cauchy(0, 1);
    R[1] ~ normal(prior_R_a, prior_R_b);

    for (i in 2 : K)
      R[i] ~ normal(R[i - 1], sigma);
  } else if (prior_R_choice == 2) {
    for(i in 1:K)
      R[i] ~ gamma_mean_sd(prior_R_a, prior_R_b);
  }

  // delay model
  // likelihood
  for (t in 1 : (N_delay - 1)) {
      if (time_onset[t] == time_onset[t + 1]) {
        int cases_observed = cases_reported[t + 1] - cases_reported[t];
        int cases_remaining = cases_true[t] - cases_reported[t];
        row_vector[2] theta_t = theta[reporting_window[t]];
        real p = p_detect(time_onset[t], time_reported[t] + 0.5,
                          time_reported[t + 1] + 0.5, theta_t,
                          is_gamma_reporting_delay);
        cases_observed ~ binomial(cases_remaining, p);
      }
  }

  // priors
  for (i in 1 : n_reporting_window) {
    theta[i, 1] ~ gamma_mean_sd(prior_theta_1_a, prior_theta_1_b);
    theta[i, 2] ~ gamma_mean_sd(prior_theta_2_a, prior_theta_2_b);
  }
}


