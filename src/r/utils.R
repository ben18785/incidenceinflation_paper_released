dgamma_mean_sd <- function(x, mu, sigma, ...) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::dgamma(x, shape, rate, ...)
}

make_generation_time_distribution <- function(generation_distribution_params) {
  w_max_days <- generation_distribution_params$w_max_days
  shape <- generation_distribution_params$shape
  rate <- generation_distribution_params$rate
  int_w <- pgamma(seq(0, w_max_days, 1), shape, rate)
  w <- diff(int_w)
  w <- w / sum(w)
  w
}

pgamma_mean_sd <- function(x, mu, sigma) {
  shape <- mu^2 / sigma^2
  rate <- mu / sigma^2
  stats::pgamma(x, shape, rate)
}

gamma_discrete_pmf <- function(day, serial_parameters){
  delay_mean <- serial_parameters$mean
  delay_sd <- serial_parameters$sd
  pgamma_mean_sd(day + 0.5, delay_mean, delay_sd) -
    pgamma_mean_sd(day - 0.5, delay_mean, delay_sd)
}

weights_series <- function(t_max, serial_parameters) {
  day_series <- seq(1, t_max, 1)
  w <- purrr::map_dbl(day_series, ~gamma_discrete_pmf(., serial_parameters))
  # due to truncation, normalise series
  w <- w / sum(w)
  w
}

quantiles_r <- function(fit) {
  R_draws <- rstan::extract(fit, "R")[[1]]
  lower <- apply(R_draws, 2, function(x) quantile(x, 0.025))
  middle <- apply(R_draws, 2, function(x) quantile(x, 0.5))
  upper <- apply(R_draws, 2, function(x) quantile(x, 0.975))
  tibble(
    lower=lower,
    middle=middle,
    upper=upper,
    Rt_index=seq_along(lower))
}


summarise_r <- function(fits) {
  fit_full <- fits$full
  fit_snapshot <- fits$snapshot
  R_full <- quantiles_r(fit_full) %>% 
    mutate(type="full")
  R_snapshot <- quantiles_r(fit_snapshot) %>% 
    mutate(type="snapshot")
  R_full %>% 
    bind_rows(R_snapshot)
}

quantiles_r_other <- function(results) {
  R_draws <- results$Rt
  max_iteration <- max(R_draws$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  R_sums <- R_draws %>% 
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(Rt_index) %>% 
    summarise(
      lower=quantile(Rt, 0.025),
      middle=quantile(Rt, 0.5),
      upper=quantile(Rt, 0.975)
    ) %>% 
    mutate(type="snapshot")
  R_sums
}

fit_both_stan_dengue <- function(dengue_processed_simple, dengue_processed,
                                 dengue_generation_ws, dengue_rt_index,
                                 end_time, niterations_dengue, model_poisson) {
  fit_both_stan(
    model_poisson,
    dengue_processed_simple %>% 
      filter(time_onset <= end_time),
    dengue_processed %>% 
      filter(time_reported <= end_time) %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported)),
    40,
    dengue_generation_ws,
    dengue_rt_index$Rt_index[1:end_time],
    I_init=rep(1, 40),
    n_iterations=niterations_dengue)
}

fit_measles <- function(time_end,
                        measles_processed,
                        measles_rt_index,
                        measles_reporting_index,
                        serial_parameters_measles,
                        niterations_measles,
                        is_negative_binomial=FALSE
) {
  fit_model(measles_processed %>% 
              filter(time_reported <= time_end) %>% 
              thin_series() %>% 
              left_join(measles_rt_index, by="time_onset") %>% 
              left_join(measles_reporting_index, by="time_onset"),
            serial_parameters_measles,
            niterations=niterations_measles,
            is_negative_binomial=is_negative_binomial)
}

fit_ebola <- function(time_remove,
                      data_sim_resurgence_ebola,
                      serial_parameters_ebola,
                      niterations_resurgence) {
  
  fit_model(data_sim_resurgence_ebola %>% 
              select(-cases_true) %>%
              filter(time_reported <= (125 - time_remove)),
            serial_parameters_ebola,
            niterations=niterations_resurgence,
            is_negative_binomial=TRUE)
}

fit_ebola_real <- function(time_end,
                      data_sim_resurgence_ebola,
                      serial_parameters_ebola,
                      niterations_resurgence) {
  
  fit_model(data_sim_resurgence_ebola %>% 
              select(-c(cases_true, date_onset)) %>%
              filter(time_reported <= time_end),
            serial_parameters_ebola,
            niterations=niterations_resurgence,
            is_negative_binomial=FALSE)
}

fit_both_stan_simple <- function(processed_simple, processed,
                                 generation_ws, rt_index,
                                 end_time, niterations, model) {
  fit_both_stan(
    model,
    processed_simple %>% 
      filter(time_onset <= end_time),
    processed %>% 
      filter(time_reported <= end_time) %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported)),
    40,
    generation_ws,
    rt_index$Rt_index[1:end_time],
    I_init=rep(1, 40),
    n_iterations=niterations)
}

fit_measles_uncertain_cases_stan <- function(time_end,
                        measles_processed,
                        measles_rt_index,
                        measles_reporting_index,
                        serial_parameters_measles,
                        niterations_measles,
                        is_negative_binomial=FALSE,
                        is_maximise=FALSE,
                        initial_cases_true=NULL,
                        initial_reporting_parameters=NULL,
                        initial_Rt=NULL,
                        initial_overdispersion=10,
                        use_stan_sampling=FALSE,
                        n_stan_warmup=n_stan_warmup,
                        n_stan_iterations=n_stan_iterations,
                        step_size=step_size
) {
  fit_model_uncertain_cases_stan(measles_processed %>% 
              filter(time_reported <= time_end) %>% 
              thin_series() %>% 
              left_join(measles_rt_index, by="time_onset") %>% 
              left_join(measles_reporting_index, by="time_onset"),
            serial_parameters_measles,
            niterations=niterations_measles,
            is_negative_binomial=is_negative_binomial,
            is_maximise=is_maximise,
            initial_cases_true=initial_cases_true,
            initial_reporting_parameters=initial_reporting_parameters,
            initial_Rt=initial_Rt,
            initial_overdispersion=initial_overdispersion,
            use_stan_sampling=use_stan_sampling,
            n_stan_warmup=n_stan_warmup,
            n_stan_iterations=n_stan_iterations,
            step_size=step_size)
}


fit_ebola_uncertain_cases_stan <- function(time_end,
                                           data_sim_resurgence_ebola,
                                           serial_parameters_ebola,
                                           niterations_resurgence,
                                           stan_model,
                                           is_negative_binomial=FALSE,
                                           is_maximise=FALSE,
                                           initial_cases_true=NULL,
                                           initial_reporting_parameters=NULL,
                                           initial_Rt=NULL,
                                           initial_overdispersion=10) {
  
  fit_model_uncertain_cases_stan(data_sim_resurgence_ebola %>% 
                                   select(-cases_true) %>%
                                   filter(time_reported <= time_end),
                                 serial_parameters_ebola,
                                 stan_model,
                                 niterations=niterations_resurgence,
                                 is_negative_binomial=is_negative_binomial,
                                 is_maximise=is_maximise,
                                 initial_cases_true=initial_cases_true,
                                 initial_reporting_parameters=initial_reporting_parameters,
                                 initial_Rt=initial_Rt,
                                 initial_overdispersion=initial_overdispersion)
}