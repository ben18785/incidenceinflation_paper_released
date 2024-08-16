
forward_simulate <- function(t_periods, current_Rt, case_history, w,
                             kappa=1000) {
  wmax <- length(w)
  I <- vector()
  for(i in 1:t_periods) {
    mu <- current_Rt * sum(w * case_history)
    I_temp <- rnbinom(1, size=kappa, mu=mu)
    I <- c(I_temp, I)
    case_history <- c(I_temp, case_history)
    case_history <- case_history[1:wmax]
  }
  
  I
}

prob_outbreak_over <- function(t_horizon, current_Rt, case_history, w,
                               kappa=1000, n_iterates=1000) {
  is_over_vector <- vector(length = n_iterates)
  for(i in seq_along(is_over_vector)) {
    I <- forward_simulate(t_horizon, current_Rt, case_history, w,
                          kappa)
    is_over_vector[i] <- if_else(I[1] == 0, 1, 0)
  }
  mean(is_over_vector)
}

get_case_history <- function(case_short_df, wmax) {
  case_history <- rev(case_short_df$cases_true)
  case_history[1:wmax]
}

prob_outbreak_over_incidence_inflation <- function(t_horizon, fit_object, w,
                                                   kappa=1000, n_iterates=1000) {
  rt_df <- fit_object$Rt
  cases_df <- fit_object$cases
  
  max_iteration <- max(rt_df$iteration)
  rt_df <- rt_df %>% 
    filter(iteration >= max_iteration / 2) %>% 
    mutate(combined=paste0(iteration, chain)) %>% 
    mutate(counter=as.numeric(as.factor(combined)))
  cases_df <- cases_df %>% 
    filter(iteration >= max_iteration / 2) %>% 
    mutate(combined=paste0(iteration, chain)) %>% 
    mutate(counter=as.numeric(as.factor(combined)))
  
  wmax <- length(w)
  rt_last_index <- max(rt_df$Rt_index)
  
  is_over_vector <- vector(length = n_iterates)
  for(i in seq_along(is_over_vector)) {
    
    id <- sample(1:max(rt_df$counter), 1)
    
    cases_short_df <- cases_df %>% filter(counter==id)
    case_history <- get_case_history(cases_short_df, wmax)
    rt_current <- rt_df %>%
      filter(counter==id) %>%
      filter(Rt_index==rt_last_index) %>% 
      pull(Rt)
    I <- forward_simulate(t_horizon, rt_current, case_history, w,
                          kappa)
    is_over_vector[i] <- if_else(I[1] == 0, 1, 0)
  }
  mean(is_over_vector)
}

forward_simulations_many <- function(
    t_horizon, current_Rt, case_history, w,
    current_time_onset,
    kappa=1000, n_iterates=100) {
  
  for(i in 1:n_iterates) {
    I <- forward_simulate(t_horizon, current_Rt, case_history, w,
                          kappa)
    df_tmp <- tibble(I=I,
                     time_onset=seq(current_time_onset + 1, current_time_onset + length(I)),
                     iteration=i)
    
    if(i == 1)
      big_df <- df_tmp
    else
      big_df <- big_df %>% bind_rows(df_tmp)
  }
  big_df
}

forward_simulations_many_incidenceinflation <- function(
    t_horizon, fit_object, w,
    current_time_onset,
    kappa=1000, n_iterates=100) {
  
  rt_df <- fit_object$Rt
  cases_df <- fit_object$cases
  
  max_iteration <- max(rt_df$iteration)
  rt_df <- rt_df %>% 
    filter(iteration >= max_iteration / 2) %>% 
    mutate(combined=paste0(iteration, chain)) %>% 
    mutate(counter=as.numeric(as.factor(combined)))
  cases_df <- cases_df %>% 
    filter(iteration >= max_iteration / 2) %>% 
    mutate(combined=paste0(iteration, chain)) %>% 
    mutate(counter=as.numeric(as.factor(combined)))
  
  wmax <- length(w)
  rt_last_index <- max(rt_df$Rt_index)
  
  for(i in 1:n_iterates) {
    
    id <- sample(1:max(rt_df$counter), 1)
    
    cases_short_df <- cases_df %>% filter(counter==id)
    case_history <- get_case_history(cases_short_df, wmax)
    rt_current <- rt_df %>%
      filter(counter==id) %>%
      filter(Rt_index==rt_last_index) %>% 
      pull(Rt)
    
    I <- forward_simulate(t_horizon, rt_current, case_history, w, kappa)
    df_tmp <- tibble(I=I,
                     time_onset=seq(current_time_onset + 1, current_time_onset + length(I)),
                     iteration=i)
    
    if(i == 1)
      big_df <- df_tmp
    else
      big_df <- big_df %>% bind_rows(df_tmp)
  }
  big_df
}
