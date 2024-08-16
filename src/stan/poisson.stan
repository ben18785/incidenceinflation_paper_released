
data {
  int N; // number of data points
  int K; // number of R segments
  int window[N]; // assigns a time point to a given segment
  int C[N]; // case series
  int wmax; // max day of generation time distribution
  row_vector[wmax] w; // generation distribution
  vector[wmax] I_history; // note the order for this should be most recent cases first
}

parameters {
  real<lower=0> R[K];
}

transformed parameters {
  
  real E_cases[N];
  real sum_w[N];
  
  {
    vector[wmax] I_temp;
    for(t in 1:N) {
      if(t == 1) {
        I_temp = I_history;
      } else if(t <= wmax) {
        int kk = wmax - t + 1;
        for(i in 1:(t - 1))
          I_temp[i] = C[t - i]; // needs to be lagged by one time point
        for(i in 1:kk)
          I_temp[i + t - 1] = I_history[i];
      } else {
        for(i in 1:wmax)
          I_temp[i] = C[t - i]; // needs to be lagged by one time point
      }
      sum_w[t] = w * I_temp;
      E_cases[t] = R[window[t]] * sum_w[t];
    }
  }
}

model {
  for(t in 1:N)
    C[t] ~ poisson(E_cases[t]);
  R ~ gamma(1.0, 0.2);
}

generated quantities {
  vector[N] log_likelihood;
  for(t in 1:N)
    log_likelihood[t] = poisson_lpmf(C[t]|E_cases[t]);
}
