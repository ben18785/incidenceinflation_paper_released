functions {
  real prob_detection_gamma(
    real time_onset,
    real time_report_start_interval,
    real time_report_end_interval,
    real alpha,
    real beta) {
    
    real time_1 = time_report_start_interval - time_onset;
    real time_2 = time_report_end_interval - time_onset;
    real survival_1 = 1 - gamma_cdf(time_1, alpha, beta);
    real survival_2 = 1 - gamma_cdf(time_2, alpha, beta);
    
    return (survival_1 - survival_2) / survival_1;
  }
  
  int bin_search(real x, int min_val, int max_val) {
    if (min_val > x || max_val < x) 
      reject("require min < x < max, found min = ", min_val, "; max = ",
             max_val, "; x = ", x);
    real y = round(x);
    int range = max_val - min_val;
    int mid_pt = min_val;
    while (1) {
      if (range == 0) 
        return mid_pt;
      range = (range + 1) / 2;
      mid_pt += y > mid_pt ? range : -range;
    }
    return min_val; // never reached
  }

}

data {
  
  int N; // number of data points
  int time_onset[N];
  real time_report_start_interval[N];
  real time_report_end_interval[N];
  int n_detected[N];
  int cumulative_detected[N];
  int time_onset_max;
  
  
  int K; // number of R segments
  int window[N]; // assigns a time point to a given segment
  int wmax; // max day of generation time distribution
  row_vector[wmax] w; // generation distribution
  vector[wmax] I_history; // note the order for this should be most recent cases first
}

parameters {
  real<lower=0> R[K];
  real<lower=0> I[time_onset_max];
  real<lower=0> alpha;
  real<lower=0> beta;
}

transformed parameters {
  
  real E_cases[time_onset_max];
  
  {
    vector[wmax] I_temp;
    for(t in 1:time_onset_max) {
      if(t == 1) {
        I_temp = I_history;
      } else if(t < wmax) {
        int kk = wmax - t;
        for(i in 1:t)
          I_temp[i] = I[t - i + 1];
        for(i in 1:kk)
          I_temp[i + t] = I_history[i];
      } else {
        for(i in 1:wmax)
          I_temp[i] = I[t - i + 1];
      }
      E_cases[t] = R[window[t]] * w * I_temp;
    }
  }
}

model {
  for(t in 1:time_onset_max) {
    int I_round[t] = bin_search(I[t], 0, 1000);
    I_round[t] ~ poisson(E_cases[t]);
  }
  for(i in 1:N) {
    int n_remaining = I_round[time_onset[i]] - cumulative_detected[i];
    real prob_detect = prob_detection_gamma(
      time_onset[i],
      time_report_start_interval[i],
      time_report_end_interval[i],
      alpha,
      beta);
    n_detected[i] ~ binomial(n_remaining, prob_detect);
  }
    
  R ~ gamma(1.0, 0.2);
}
