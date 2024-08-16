
plot_cases_vs_true <- function(results, df_true) {
  
  cases_df <- results$cases
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true)
  
  cases_df %>% 
    ggplot(aes(x=time_onset, y=cases_true)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(colour="black") +
    xlab("Onset date") +
    ylab("Cases")
}

plot_cases_vs_true_series <- function(cases_df_combined, df_true, include_reporting=FALSE, transform_type="log") {
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true)
  
  g <- cases_df %>% 
    ggplot(aes(x=time_onset)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(aes(y=cases_true), colour="black") +
    xlab("Onset date") +
    ylab("Cases") +
    facet_wrap(~max_time)
  
  if(include_reporting)
    g <- g + geom_line(aes(y=cases_reported), colour="#FF8C00")
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt(limits=c(0, 100))
  g
}

plot_cases_vs_true_series_2 <- function(cases_df_combined, df_true, df_full, include_reporting=FALSE, transform_type="log",
                                        ylimits=c(0, 52), min_date=NULL) {
  
  lookup <- df_full %>% 
    select(date_onset, time_onset) %>% 
    unique()
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true) %>% 
    left_join(lookup)
  
  final_date <- max(cases_df$date_onset)
  max_times <- as.numeric(as.character(sort(unique(cases_df$max_time))))
  dates <- vector(length = length(max_times))
  for(i in seq_along(max_times)) {
    dates[i] <- as.character(final_date - max_times[i])
  }
  df_dates_lookup <- tibble(max_time=max_times, observation_date=dates) %>% 
    mutate(observation_date=as.Date(observation_date)) %>% 
    mutate(max_time=as.character(max_time)) %>% 
    mutate(observation_date=format(observation_date, "%d %b %Y"))
  cases_df <- cases_df %>% 
    left_join(df_dates_lookup) %>% 
    mutate(observation_date=as.factor(observation_date)) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    mutate(observation_date=fct_rev(fct_reorder(observation_date, max_time)))
    
  if(!is.null(min_date))
    cases_df <- cases_df %>% filter(date_onset>=min_date)
  
  
  if(include_reporting) {
    
    g <- cases_df %>% 
      ggplot(aes(x=date_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(aes(y=cases_reported), colour="#FF8C00") +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset date") +
      ylab("Cases") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      facet_wrap(~observation_date)
      
    } else {
      g <- cases_df %>% 
        ggplot(aes(x=date_onset)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
        geom_line(aes(y=middle), colour="#1B9E77") +
        geom_line(aes(y=cases_true), colour="black") +
        xlab("Onset date") +
        ylab("Cases") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
        facet_wrap(~observation_date)
    }
      
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt(limits=ylimits)
  g
}

plot_cases_vs_true_series_3 <- function(cases_df_combined, df_true, df_full, observation_dates, include_reporting=FALSE, transform_type="log",
                                        ylimits=c(0, 52), min_date=NULL) {
  
  lookup <- df_full %>% 
    select(date_onset, time_onset) %>% 
    unique()
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true) %>% 
    left_join(lookup)
  
  final_date <- max(cases_df$date_onset)
  df_dates_lookup <- observation_dates %>% 
    mutate(observation_date=as.Date(observation_point)) %>% 
    mutate(max_time=as.character(max_time)) %>% 
    mutate(observation_date=format(observation_date, "%d %b %Y"))
  cases_df <- cases_df %>% 
    left_join(df_dates_lookup) %>% 
    mutate(observation_date=as.factor(observation_date)) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    mutate(observation_date=fct_reorder(observation_date, max_time))
  
  if(!is.null(min_date))
    cases_df <- cases_df %>% filter(date_onset>=min_date)
  
  
  if(include_reporting) {
    
    g <- cases_df %>% 
      ggplot(aes(x=date_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(aes(y=cases_reported), colour="#FF8C00") +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset date") +
      ylab("Cases") +
      scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
      facet_wrap(~observation_date)
    
  } else {
    g <- cases_df %>% 
      ggplot(aes(x=date_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset date") +
      ylab("Cases") +
      scale_x_date(date_breaks = "2 month2", date_labels = "%b %Y") +
      facet_wrap(~observation_date)
  }
  
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt(limits=ylimits)
  g
}

plot_cases_vs_true_series_one_line <- function(cases_df_combined, df_true, df_full) {
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  lookup <- df_full %>% 
    select(time_onset, onset_week) %>% 
    unique()
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true) %>% 
    left_join(lookup) %>% 
    group_by(max_time) %>% 
    mutate(observation_point=max(onset_week)) %>% 
    mutate(observation_point=format(observation_point, "%d %b %Y")) %>% 
    ungroup() %>% 
    mutate(observation_point=as.factor(observation_point)) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    mutate(observation_point=fct_rev(fct_reorder(observation_point, max_time)))
  
  g <- cases_df %>% 
    filter(time_onset >= 60) %>% 
    ggplot(aes(x=onset_week)) +
    geom_line(aes(y=cases_reported, group=as.factor(observation_point), colour=as.factor(observation_point)),
              linetype=2) +
    geom_ribbon(aes(ymin=lower, ymax=upper, group=as.factor(observation_point), fill=as.factor(observation_point)), alpha=0.4) +
    geom_line(aes(y=middle, group=as.factor(observation_point), colour=as.factor(observation_point))) +
    geom_line(aes(y=cases_true), colour="black") +
    scale_color_brewer("Observation\ndate", palette = "Set1") +
    scale_fill_brewer("Observation\ndate", palette = "Set1") +
    xlab("Onset date") +
    ylab("Cases") +
    scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
    theme(
      text=element_text(size=16)
    )
  
  g
}

plot_cases_vs_true_series_1 <- function(cases_df_combined, df_true, include_reporting=FALSE, transform_type="log") {
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  max_time_lookup <- cases_df_combined %>%
    select(-c(cases_true, iteration, chain)) %>% 
    unique()
  reporting_df <- df_true %>%
    left_join(max_time_lookup) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    filter(time_reported <= max_time) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      cases_reported=last(cases_reported)
    )
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true)
  
  if(include_reporting)
    g <- cases_df %>% 
    ggplot(aes(x=time_onset)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(data=reporting_df,
              aes(y=cases_reported), colour="#FF8C00") +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(aes(y=cases_true), colour="black") +
    xlab("Onset time") +
    ylab("Cases") +
    facet_wrap(~max_time)
  else
    g <- cases_df %>% 
    ggplot(aes(x=time_onset)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(aes(y=cases_true), colour="black") +
    xlab("Onset time") +
    ylab("Cases") +
    facet_wrap(~max_time)
    
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt()
  g
}

plot_rt_vs_true <- function(results, rt_true, df_true, threshold_onset_time=1) {
  
  res_df <- results$Rt
  max_iteration <- max(res_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  res_df <- res_df %>% 
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(Rt_index) %>% 
    summarise(
      lower=quantile(Rt, 0.025),
      middle=quantile(Rt, 0.5),
      upper=quantile(Rt, 0.975)
    )
  
  full_df <- tibble(
    time_onset=seq_along(rt_true),
    Rt=rt_true
  ) %>% 
    left_join(df_true, by="time_onset") %>% 
    left_join(res_df, by="Rt_index")
  
  full_df %>% 
    filter(time_onset >= threshold_onset_time) %>% 
    ggplot(aes(x=time_onset, y=Rt)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(colour="black") +
    xlab("Onset date") +
    ylab(TeX("$R_t$")) +
    geom_hline(yintercept = 1, linetype=2)
}

plot_rt_vs_true_series <- function(rt_combined, rt_true, df_true, threshold_onset_time=1) {
  
  res_df <- rt_combined
  max_iteration <- max(res_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  res_df <- res_df %>% 
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, Rt_index) %>% 
    summarise(
      lower=quantile(Rt, 0.025),
      middle=quantile(Rt, 0.5),
      upper=quantile(Rt, 0.975)
    )
  
  full_df <- tibble(
    time_onset=seq_along(rt_true),
    Rt=rt_true
  ) %>% 
    left_join(df_true, by="time_onset") %>% 
    left_join(res_df, by="Rt_index")
  
  full_df %>% 
    filter(time_onset >= threshold_onset_time) %>% 
    ggplot(aes(x=time_onset, y=Rt)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    geom_line(colour="black") +
    xlab("Onset date") +
    ylab(TeX("$R_t$")) +
    geom_hline(yintercept = 1, linetype=2) +
    facet_wrap(~max_time)
}


plot_delay_true <- function(results, reporting_parameters_true) {
  
  mu <- reporting_parameters_true$mean
  sd <- reporting_parameters_true$sd
  x <- seq(0, 20, length.out=100)
  pdf_true <- dgamma_mean_sd(x, mu, sd)
  df_true <- tibble(
    x=x, pdf=pdf_true
  )
  
  rep_df <- results$reporting
  max_iteration <- max(rep_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  rep_df <- rep_df %>% 
    filter(iteration >= post_warmup_iteration) 
  
  ndraws <- 100
  
  for(i in 1:ndraws) {
    index <- sample(nrow(rep_df), 1)
    row_draws <- rep_df[index, ]
    pdf_draw <- dgamma_mean_sd(x, row_draws$mean, row_draws$sd)
    df_draw <- tibble(
      x=x, pdf=pdf_draw, draw=i
    )
    if(i == 1)
      df_draws <- df_draw
    else
      df_draws <- df_draws %>% bind_rows(df_draw)
  }
  df_draws <- df_draws %>% 
    group_by(x) %>% 
    summarise(
      lower=quantile(pdf, 0.01),
      middle=quantile(pdf, 0.5),
      upper=quantile(pdf, 0.99)
    )
  
  ggplot(df_true, aes(x=x)) +
    geom_ribbon(data=df_draws,
                aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=pdf)) +
    xlab("Reporting delay, days") +
    ylab("Probability density")
}

construct_pdf <- function(mean, sd, x) {
  
  pdf_true <- dgamma_mean_sd(x, mean, sd)
  df_true <- tibble(
    x=x, pdf=pdf_true
  )
  df_true
}

sample_summarise_densities <- function(rep_df, ndraws, x) {
  
  for(i in 1:ndraws) {
    index <- sample(nrow(rep_df), 1)
    row_draws <- rep_df[index, ]
    pdf_draw <- dgamma_mean_sd(x, row_draws$mean, row_draws$sd)
    df_draw <- tibble(
      x=x, pdf=pdf_draw, draw=i
    )
    if(i == 1)
      df_draws <- df_draw
    else
      df_draws <- df_draws %>% bind_rows(df_draw)
  }
  
  df_draws %>% 
    group_by(x) %>% 
    summarise(
      lower=quantile(pdf, 0.01),
      middle=quantile(pdf, 0.5),
      upper=quantile(pdf, 0.99)
    )
}


plot_delay_true_series <- function(reporting_delay_combined, reporting_parameters_true) {
  
  x <- seq(0, 20, length.out=100)
  pdf_true <- construct_pdf(reporting_parameters_true$mean,
                            reporting_parameters_true$sd,
                            x)
  
  max_iteration <- max(reporting_delay_combined$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  rep_df <- reporting_delay_combined %>% 
    filter(iteration >= post_warmup_iteration) 
  
  max_times <- unique(reporting_delay_combined$max_time)
  means_at_median_posteriors <- vector(length = length(max_times))
  for(i in seq_along(max_times)) {
    rep_df_tmp <- rep_df %>% filter(max_time == max_times[i])
    tmp <- sample_summarise_densities(rep_df_tmp, 100, x) %>% 
      mutate(max_time=max_times[i])
    if(i == 1)
      df_sim <- tmp
    else
      df_sim <- df_sim %>% bind_rows(tmp)
    # means_at_median_posteriors[i] <- 
  }
  
  ggplot(pdf_true, aes(x=x)) +
    geom_ribbon(data=df_sim,
                aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=pdf)) +
    xlab("Reporting delay, days") +
    ylab("Probability density") +
    facet_wrap(~max_time)
}

# works only for one max_time
plot_delay_true_series_simple <- function(reporting_delay_combined, reporting_parameters_true) {
  
  x <- seq(0, 20, length.out=100)
  pdf_true <- construct_pdf(reporting_parameters_true$mean,
                            reporting_parameters_true$sd,
                            x)
  
  max_iteration <- max(reporting_delay_combined$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  rep_df <- reporting_delay_combined %>% 
    filter(iteration >= post_warmup_iteration) 
  
  max_times <- unique(reporting_delay_combined$max_time)
  for(i in seq_along(max_times)) {
    rep_df_tmp <- rep_df %>% filter(max_time == max_times[i])
    tmp <- sample_summarise_densities(rep_df_tmp, 100, x) %>% 
      mutate(max_time=max_times[i])
    if(i == 1)
      df_sim <- tmp
    else
      df_sim <- df_sim %>% bind_rows(tmp)
  }
  
  ggplot(pdf_true, aes(x=x)) +
    geom_ribbon(data=df_sim,
                aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=pdf)) +
    xlab("Reporting delay, days") +
    ylab("Probability density")
}

plot_delay_time_varying <- function(df_delays, x=seq(0, 40, length.out=100)) {
  
  roman_num <- c("(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)", "(vii)", "(viii)", "(ix)", "(x)",
                  "(xi)", "(xii)", "(xiii)", "(xiv)", "(xv)", "(xvi)", "(xvii)", "(xviii)", "(xix)", "(xx)")
  
  toLowerRoman <- function(x) {
    if (x < 1 | x > 20) {
      stop("Input must be between 1 and 20")
    }
    return(roman_num[x])
  }
  
  max_iteration <- max(df_delays$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  rep_df <- df_delays %>% 
    filter(iteration >= post_warmup_iteration) 
  
  rep_indices <- unique(df_delays$reporting_piece_index)
  for(i in seq_along(rep_indices)) {
    rep_df_tmp <- rep_df %>%
      filter(reporting_piece_index == rep_indices[i])
    tmp <- sample_summarise_densities(rep_df_tmp, 100, x) %>% 
      mutate(reporting_piece_index=rep_indices[i]) %>% 
      mutate(label=toLowerRoman(reporting_piece_index)) %>% 
      mutate(label=as.factor(label))
    if(i == 1)
      df_sim <- tmp
    else
      df_sim <- df_sim %>% bind_rows(tmp)
  }
  
  
  ggplot(df_sim, aes(x=x)) +
    geom_ribbon(data=df_sim,
                aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
    geom_line(aes(y=middle), colour="#1B9E77") +
    xlab("Reporting delay, days") +
    ylab("Probability density") +
    facet_wrap(~label)
}

plot_cases_simple <- function(
    cases_df_combined, df_true, include_reporting=TRUE, transform_type="none") {
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  max_time_lookup <- cases_df_combined %>%
    select(-c(cases_true, iteration, chain)) %>% 
    unique()
  reporting_df <- df_true %>%
    left_join(max_time_lookup) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    filter(time_reported <= max_time) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      cases_reported=last(cases_reported)
    )
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true)
  
  
  if(include_reporting) {
    g <- cases_df %>% 
      filter(max_time=="100") %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(data=reporting_df %>% filter(max_time=="100"),
                aes(y=cases_reported), colour="#FF8C00") +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset date") +
      ylab("Cases")
  } else{
    g <- cases_df %>% 
      filter(max_time=="100") %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset date") +
      ylab("Cases")
  }
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt(limits=c(0, 100))
  g
}

plot_cases_simple_schematic <- function(
    cases_df_combined, df_true, include_reporting=TRUE, transform_type="none") {
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  max_time_lookup <- cases_df_combined %>%
    select(-c(cases_true, iteration, chain)) %>% 
    unique()
  reporting_df <- df_true %>%
    left_join(max_time_lookup) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    filter(time_reported <= max_time) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      cases_reported=last(cases_reported)
    ) %>% 
    filter(time_onset>=25) %>% 
    mutate(time_onset=time_onset-25+1)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true) %>% 
    ungroup() %>% 
    filter(time_onset>=25) %>% 
    mutate(time_onset=time_onset-25+1)
  
  
  if(include_reporting) {
    g <- cases_df %>% 
      filter(max_time=="100") %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(data=reporting_df %>% filter(max_time=="100"),
                aes(y=cases_reported), colour="#FF8C00") +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset day") +
      ylab("Cases")
  } else{
    g <- cases_df %>% 
      filter(max_time=="100") %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_true), colour="black") +
      xlab("Onset day") +
      ylab("Cases")
  }
  
  if(transform_type=="log")
    g <- g + scale_y_log10()
  else if(transform_type=="sqrt")
    g <- g + scale_y_sqrt(limits=c(0, 100))
  g
}

process_cases_data <- function(cases_df_combined, df_true, df_full, observation_dates) {
  
  lookup <- df_full %>% 
    select(date_onset, time_onset) %>% 
    unique()
  
  cases_df <- cases_df_combined
  max_iteration <- max(cases_df$iteration)
  post_warmup_iteration <- floor(max_iteration / 2)
  
  cases_df <- cases_df %>%
    filter(iteration >= post_warmup_iteration) %>% 
    group_by(max_time, time_onset) %>% 
    summarise(
      lower=quantile(cases_true, 0.025),
      middle=quantile(cases_true, 0.5),
      upper=quantile(cases_true, 0.975)
    ) %>% 
    left_join(df_true) %>% 
    left_join(lookup)
  
  final_date <- max(cases_df$date_onset)
  df_dates_lookup <- observation_dates %>% 
    mutate(observation_date=as.Date(observation_point)) %>% 
    mutate(max_time=as.character(max_time)) %>% 
    mutate(observation_date=format(observation_date, "%d %b %Y"))
  cases_df <- cases_df %>% 
    left_join(df_dates_lookup) %>% 
    mutate(observation_date=as.factor(observation_date)) %>% 
    mutate(max_time=as.numeric(as.character(max_time))) %>% 
    mutate(observation_date=fct_reorder(observation_date, max_time))
  
  cases_df
}