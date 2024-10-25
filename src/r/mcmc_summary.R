
mcmc_summary <- function(mcmc_results) {
  
  mcmc_results1 <- incidenceinflation::convert_results_to_posterior_format(mcmc_results)
  cases_df <- mcmc_results1$cases
  Rt_df <- mcmc_results1$Rt
  rep_df <- mcmc_results1$reporting
  
  is_negative_binomial <- FALSE
  if("overdispersion" %in% names(mcmc_results1)) {
    is_negative_binomial <- TRUE
    overdispersion_df <- mcmc_results1$overdispersion
  }
  
  combined_df <- rep_df %>% 
    left_join(Rt_df) %>% 
    left_join(cases_df)
  
  if(is_negative_binomial)
    combined_df <- combined_df %>% 
      left_join(overdispersion_df)
  
  posterior::summarise_draws(combined_df)
}

mcmc_summary_stan_uncertain_cases <- function(mcmc_results) {
  
  mcmc_results1 <- incidenceinflationstan::convert_results_to_posterior_format(mcmc_results)
  cases_df <- mcmc_results1$cases
  Rt_df <- mcmc_results1$Rt
  rep_df <- mcmc_results1$reporting
  
  is_negative_binomial <- FALSE
  if("other" %in% names(mcmc_results1)) {
    is_negative_binomial <- TRUE
    overdispersion_df <- mcmc_results1$other
  }
  
  combined_df <- rep_df %>% 
    left_join(Rt_df) %>% 
    left_join(cases_df)
  
  if(is_negative_binomial)
    combined_df <- combined_df %>% 
    left_join(overdispersion_df)
  
  posterior::summarise_draws(combined_df)
}