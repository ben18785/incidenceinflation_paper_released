
data {
  int N; // number of data points
  int K; // number of R segments
  int window[N]; // assigns a time point to a given segment
  int C[N]; // case series
  int wmax; // max day of generation time distribution
  row_vector[wmax] w; // generation distribution
}

parameters {
  real<lower=0> R[K];
  real<lower=0> sigma;
}

transformed parameters {
  
  real E_cases[N - 1];
  
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
  for(t in 2:N)
    C[t] ~ poisson(E_cases[t - 1]);
  R[1] ~ gamma(1.0, 0.2);
  for(i in 2:K) {
    R[i] ~ normal(R[i - 1], sigma);
  }
  sigma ~ normal(0, 1);
}
