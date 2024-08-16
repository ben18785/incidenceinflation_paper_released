undetected_prob <- function(t_1, t_onset, reporting_parameters){
  if(t_1 < t_onset)
    stop("t_1 must equal or exceed t_onset")
  delay_mean <- reporting_parameters$mean
  delay_sd <- reporting_parameters$sd
  1 - pgamma_mean_sd(t_1 - t_onset, delay_mean, delay_sd)
}

detected_after_unobserved_prob <- function(day_2, day_1, day_onset,
                                           reporting_parameters){
  if(day_2 < day_1)
    stop("second observation date must be at same time or after first")
  (undetected_prob(day_1 + 0.5, day_onset, reporting_parameters) -
      undetected_prob(day_2 + 0.5, day_onset, reporting_parameters)) /
    undetected_prob(day_1 + 0.5, day_onset, reporting_parameters)
}

observed_cases_single <- function(cases_observed, cases_true, day_2, day_1, day_onset,
                                  reporting_parameters, kappa){
  if(cases_true < cases_observed)
    stop("true case count must exceed reported")
  I_remaining <- cases_true - cases_observed
  if(I_remaining == 0)
    cases <- 0
  else {
    p <- detected_after_unobserved_prob(day_2, day_1, day_onset,
                                        reporting_parameters)
    theta <- rbeta(1, p * kappa, (1 - p) * kappa)
    cases <- stats::rbinom(1, I_remaining, theta)
  }
  cases
}

observed_cases_trajectory <- function(cases_true, days_reporting, day_onset,
                                      reporting_parameters, kappa){
  
  I_observed <- vector(length = length(days_reporting))
  for(i in 1:length(days_reporting)){
    if(i == 1) {
      I_previous_obs <- 0
      day_previous_report <- day_onset
    } else {
      I_previous_obs <- I_observed[i - 1]
      day_previous_report <- days_reporting[i - 1]
    }
    new_cases <- observed_cases_single(I_previous_obs, cases_true,
                                       days_reporting[i],
                                       day_previous_report,
                                       day_onset,
                                       reporting_parameters,
                                       kappa)
    I_observed[i] <- I_previous_obs + new_cases
  }
  I_observed
}


create_reporting_from_single_parameters_df <- function(time_onsets,
                                                       reporting_parameters) {
  
  dplyr::tibble(time_onset=time_onsets,
                mean=reporting_parameters$mean,
                sd=reporting_parameters$sd
  )
}

observed_cases <- function(cases_true, reporting_parameters,
                           kappa,
                           days_max_follow_up=30){
  
  if(methods::is(reporting_parameters, "list"))
    reporting_parameters <- create_reporting_from_single_parameters_df(
      seq_along(cases_true),
      reporting_parameters)
  
  if(nrow(reporting_parameters) != length(cases_true))
    stop("Number of rows in reporting parameters tibble must match number of cases.")
  
  d_max <- length(cases_true)
  for(t in 1:d_max){
    
    a_max <- min(c(d_max, t + days_max_follow_up))
    obs_time <- seq(t, a_max, 1)
    reporting_parameters_current <- list(mean=reporting_parameters$mean[t],
                                         sd=reporting_parameters$sd[t])
    cases_obs_trajec <- observed_cases_trajectory(
      cases_true=cases_true[t],
      days_reporting=obs_time,
      day_onset=t,
      reporting_parameters_current,
      kappa=kappa)
    
    # stack into data frame containing identifying info
    I_obs_single_onset <- dplyr::tibble(
      time_onset=rep(t, length(cases_obs_trajec)),
      time_reported=obs_time,
      cases_reported=cases_obs_trajec,
      cases_true=cases_true[t])
    if(t == 1)
      cases_obs <- I_obs_single_onset
    else
      cases_obs <- cases_obs %>%
      dplyr::bind_rows(I_obs_single_onset)
  }
  cases_obs
}