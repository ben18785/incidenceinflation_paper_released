
days_total <- 100
v_Rt <- c(rep(1.5, 40), rep(1.2, 20), seq(1, 1.5, length.out=40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(mean=3, sd=5)
kappa <- 1000
cases <- incidenceinflation::true_cases(days_total, Rt_function, kappa, s_params,
                    initial_parameters=list(mean=30, length=30))
df <- incidenceinflation::observed_cases(cases, r_params, days_max_follow_up=40)
t <- 80
df %>% 
  filter(time_onset <= t) %>% 
  filter(time_reported >= t) %>% 
  ggplot(aes(x=time_onset, y=cases_reported)) +
  geom_line(aes(colour=time_reported,
                group=time_reported)) +
  geom_point(aes(y=cases_true), colour="black") +
  geom_line(aes(y=cases_true), colour="black") +
  scale_colour_viridis_c(begin = 0.2)


model <- rstan::stan_model("src/stan/poisson_rw.stan")

time_reporteds <- seq(t, days_total, 1)
window <- unlist(map(seq(1, 80, 1), ~rep(., 1)))

for(i in seq_along(time_reporteds)) {
  
  print(i)
  
  df_tmp <- df %>% 
    filter(time_onset <= t) %>% 
    filter(time_reported <= time_reporteds[i]) %>% 
    group_by(time_onset) %>% 
    summarise(cases_reported=max(cases_reported))
  
  data_stan <- list(
    N=nrow(df_tmp),
    K=max(window),
    C=df_tmp$cases_reported,
    window=window,
    wmax=40,
    w=weights_series(40, s_params)
  )
  
  opt <- rstan::optimizing(model, data=data_stan, as_vector=FALSE)
  df_r_tmp <- tibble(R=opt$par$R,
                     index=seq(1, data_stan$K, 1),
                     time_reported=time_reporteds[i])
  
  if(i == 1)
    df_R <- df_r_tmp
  else
    df_R <- df_R %>% bind_rows(df_r_tmp)
}

lookup <- tibble(
  index=window,
  time_onset=seq_along(index)
)

df_R %>% 
  left_join(lookup, by="index") %>% 
  filter(time_onset >= 50) %>% 
  ggplot(aes(x=time_onset, y=R, colour=time_reported, group=time_reported)) +
  geom_line() +
  scale_colour_viridis_c(begin = 0.2)


# try a different thing: vary the reporting delay and visualise cases
cases <- incidenceinflation::true_cases(days_total, Rt_function, kappa, s_params,
                                        initial_parameters=list(mean=30, length=30))

r_mus <- seq(1, 20, length.out=20)
r_params <- list(mean=3, sd=1)
for(i in seq_along(r_mus)) {
  r_params$mean <- r_mus[i]
  df <- incidenceinflation::observed_cases(cases, r_params, days_max_follow_up=40) %>% 
    mutate(delay=r_mus[i])
  tmp <- df %>% 
    group_by(time_onset) %>% 
    summarise(C=last(cases_reported)) %>% 
    ungroup() %>% 
    mutate(delay=r_mus[i])
  if(i == 1) {
    big_df <- tmp
    overall_df <- df 
  } else {
    big_df <- big_df %>% bind_rows(tmp)
    overall_df <- overall_df %>% bind_rows(df)
  }
}

t_threshold <- 25
df_true <- tibble(C=cases,
                  time_onset=seq_along(C)) %>% 
  filter(time_onset >= t_threshold) %>% 
  mutate(time_onset=seq_along(C))
df_cases_processed <- big_df %>% 
  filter(time_onset >= t_threshold) %>% 
  ungroup() %>% 
  group_by(delay) %>% 
  mutate(time_onset=seq_along(C)) %>% 
  mutate(type="Cases")
 df_cases_processed %>% 
  ggplot(aes(x=time_onset, y=C)) +
  geom_line(aes(colour=delay, group=as.factor(delay))) +
  geom_line(data=df_true) +
  scale_colour_viridis_c()

window <- unlist(map(seq(1, 100, 1), ~rep(., 1)))
lookup <- tibble(
  index=window,
  time_onset=seq_along(index)
)

for(i in seq_along(r_mus)) {
  
  print(i)
  
  df_tmp <- overall_df %>% 
    filter(delay==r_mus[i]) %>% 
    group_by(time_onset) %>% 
    summarise(cases_reported=max(cases_reported))
  
  data_stan <- list(
    N=nrow(df_tmp),
    K=max(window),
    C=df_tmp$cases_reported,
    window=window,
    wmax=40,
    w=weights_series(40, s_params)
  )
  
  opt <- rstan::optimizing(model, data=data_stan, as_vector=FALSE)
  df_r_tmp <- tibble(R=opt$par$R,
                     index=seq(1, data_stan$K, 1),
                     delay=r_mus[i])
  
  if(i == 1)
    df_R <- df_r_tmp
  else
    df_R <- df_R %>% bind_rows(df_r_tmp)
}

df_R_true <- tibble(
  R=v_Rt,
  time_onset=seq_along(R)
) %>% 
  filter(time_onset >= t_threshold) %>% 
  mutate(time_onset=seq_along(R))

df_R %>% 
  left_join(lookup, by="index") %>% 
  filter(time_onset >= t_threshold) %>% 
  group_by(delay) %>% 
  mutate(time_onset=seq_along(R)) %>% 
  ggplot(aes(x=time_onset, y=R)) +
  geom_line(aes(colour=delay, group=delay)) +
  geom_line(data=df_R_true, linetype=2) +
  scale_colour_viridis_c(begin = 0.2) +
  ylab(TeX("$R_t$"))

# Combining in a single plot
df_R_processed <- df_R %>% 
  left_join(lookup, by="index") %>% 
  filter(time_onset >= t_threshold) %>% 
  group_by(delay) %>% 
  mutate(time_onset=seq_along(R)) %>% 
  mutate(type="R[t]")

df_both <- df_cases_processed %>% 
  rename(value=C) %>% 
  bind_rows(df_R_processed %>% rename(value=R))

df_both_true <- df_true %>% 
  mutate(type="Cases") %>% 
  rename(value=C) %>% 
  bind_rows(
    df_R_true %>% 
      mutate(type="R[t]") %>% 
      rename(value=R))

# top panel
df_both %>% 
  ggplot(aes(x=time_onset, y=value)) +
  geom_line(aes(colour=delay, group=delay)) +
  geom_line(data=df_both_true, linetype=2) +
  scale_colour_viridis_c("Mean\ndelay", trans="reverse") +
  facet_wrap(~type, scales="free_y", labeller = label_parsed) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("Onset date")

# bottom panels
i <- length(r_mus)
df_tmp <- overall_df %>% 
  filter(delay==r_mus[i]) %>% 
  mutate(reporting_delay=time_reported-time_onset)

df_tmp %>%
  filter(time_onset >= t_threshold) %>% 
  ggplot(aes(x=time_reported,
             y=cases_reported,
             group=as.factor(time_onset),
             colour=reporting_delay)) +
  geom_line() +
  xlab("Onset date & report date") +
  ylab("Cases reported") +
  scale_colour_viridis_c("Delay", trans="reverse") +
  scale_y_sqrt() +
  theme(legend.position = c(0.2, 0.65))

# density panel: no need for this to be exactly to scale since it is only representative
x <- seq(0, 40, 0.1)
dens <- dgamma(x, 40, 2)
dens_df_true <- tibble(delay=x, y=dens)

nreps <- 100
for(i in 1:nreps) {
  shape <- rnorm(1, 40.5, 1)
  rate <- 2
  dens <- dgamma(x, shape, rate)
  dens_df <- tibble(delay=x, y=dens, iterate=i)
  if(i == 1)
    dens_sampled_df <- dens_df
  else
    dens_sampled_df <- dens_sampled_df %>% bind_rows(dens_df)
}

dens_sampled_df %>% 
  ggplot(aes(x=delay, y=y)) +
  geom_line(aes(group=as.factor(iterate)), alpha=0.1, colour="#1B9E77") +
  geom_line(data=dens_df_true, linetype=2) +
  xlab("Reporting delay, days") +
  ylab("Density")

