
generate_data <- function(v_Rt,
                          serial_parameters,
                          reporting_parameters,
                          Rt_window=10,
                          days_total=100,
                          kappa=1000) {
  
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  
  Rt_index_count <- floor(days_total / Rt_window)
  Rt_indices <- unlist(map(seq(1, Rt_index_count), ~rep(., Rt_window)))
  days_leftover <- days_total - length(Rt_indices)
  if(days_leftover > 0)
    Rt_indices <- c(Rt_indices, rep(Rt_index_count, days_leftover))
  stopifnot(length(Rt_indices)==days_total)
  df_Rt_index <- tibble(
    time_onset=seq(1, days_total, 1),
    Rt_index=Rt_indices
  )
  df <- incidenceinflation::generate_snapshots(days_total, Rt_function,
                           serial_parameters, reporting_parameters,
                           kappa=kappa, thinned=T) %>% 
    left_join(df_Rt_index, by="time_onset")
  df
}