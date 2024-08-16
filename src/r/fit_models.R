
fit_model <- function(df,
                      serial_parameters,
                      niterations=400,
                      nchains=4, 
                      is_negative_binomial=FALSE) {
  
  # number of reporting distributions
  num_reporting_pieces <- 1
  if("reporting_piece_index" %in% colnames(df)) {
    num_reporting_pieces <- max(df$reporting_piece_index)
  }
  
  # initial values of parameter estimates
  initial_cases_true=df %>% 
    group_by(time_onset) %>% 
    summarise(cases_true=last(cases_reported))
  
  initial_reporting_parameters = list(mean=10, sd=3)
  
  Rt_index_count <- max(df$Rt_index)
  initial_Rt <- tibble(
    Rt_index=seq(1, Rt_index_count, 1),
    Rt=1
  )
  
  # priors
  Rt_prior <- list(shape=1, rate=0.2) ## allows ~81% of probability mass with Rt>1
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=10,
                                mean_sigma=15,
                                sd_mu=5,
                                sd_sigma=15), ## uninformative reporting prior
                 max_cases=max(df$cases_reported) * 10) ## wide range of cases
  
  if(is_negative_binomial) {
    print("Fitting using a NB renewal model...")
    priors$overdispersion <- list(mean=10, sd=10)
  }
  
  if(num_reporting_pieces > 1) {
    initial_reporting_parameters <- tibble(
      reporting_piece_index=seq(1, num_reporting_pieces, 1),
      mean=initial_reporting_parameters$mean,
      sd=initial_reporting_parameters$sd
    )
  }
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  res <- incidenceinflation::mcmc(niterations=niterations,
              df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=4, is_parallel=TRUE,
              is_negative_binomial=is_negative_binomial,
              initial_overdispersion=10,
              overdispersion_metropolis_sd=1)
  stopCluster(cl)
  
  res
}

opt_model <- function(df,
                      serial_parameters,
                      niterations=400,
                      nchains=4,
                      is_negative_binomial=FALSE) {
  
  # number of reporting distributions
  num_reporting_pieces <- 1
  if("reporting_piece_index" %in% colnames(df)) {
    num_reporting_pieces <- max(df$reporting_piece_index)
  }
  
  # initial values of parameter estimates
  initial_cases_true=df %>% 
    group_by(time_onset) %>% 
    summarise(cases_true=last(cases_reported))
  
  initial_reporting_parameters = list(mean=10, sd=3)
  
  Rt_index_count <- max(df$Rt_index)
  initial_Rt <- tibble(
    Rt_index=seq(1, Rt_index_count, 1),
    Rt=1
  )
  
  # priors
  Rt_prior <- list(shape=1, rate=0.2) ## allows ~81% of probability mass with Rt>1
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=10,
                                mean_sigma=15,
                                sd_mu=5,
                                sd_sigma=15), ## uninformative reporting prior
                 max_cases=max(df$cases_reported) * 10) ## wide range of cases
  
  if(is_negative_binomial) {
    print("Fitting using a NB renewal model...")
    priors$overdispersion <- list(mean=10, sd=10)
  }
  
  if(num_reporting_pieces > 1) {
    initial_reporting_parameters <- tibble(
      reporting_piece_index=seq(1, num_reporting_pieces, 1),
      mean=initial_reporting_parameters$mean,
      sd=initial_reporting_parameters$sd
    )
  }
  
  res <- incidenceinflation::mcmc(niterations=niterations,
                                  df,
                                  priors,
                                  serial_parameters,
                                  initial_cases_true,
                                  initial_reporting_parameters,
                                  initial_Rt,
                                  reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                                  serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE,
                                  nchains=1, is_parallel=TRUE,
                                  is_negative_binomial=is_negative_binomial,
                                  initial_overdispersion=10,
                                  overdispersion_metropolis_sd=1)
  
  res
}

fit_model_inits <- function(df,
                      serial_parameters,
                      inits,
                      niterations=400,
                      nchains=4,
                      is_negative_binomial=FALSE) {
  
  # number of reporting distributions
  num_reporting_pieces <- 1
  if("reporting_piece_index" %in% colnames(df)) {
    num_reporting_pieces <- max(df$reporting_piece_index)
  }
  
  # initial values of parameter estimates
  initial_cases_true = inits$cases
  
  initial_reporting_parameters = inits$reporting
  
  Rt_index_count <- max(df$Rt_index)
  initial_Rt <- inits$Rt
  
  # priors
  Rt_prior <- list(shape=1, rate=0.2) ## allows ~81% of probability mass with Rt>1
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=10,
                                mean_sigma=15,
                                sd_mu=5,
                                sd_sigma=15), ## uninformative reporting prior
                 max_cases=max(df$cases_reported) * 10) ## wide range of cases
  
  if(is_negative_binomial) {
    print("Fitting using a NB renewal model...")
    priors$overdispersion <- list(mean=10, sd=10)
  }
  
  if(num_reporting_pieces > 1) {
    initial_reporting_parameters <- tibble(
      reporting_piece_index=seq(1, num_reporting_pieces, 1),
      mean=initial_reporting_parameters$mean,
      sd=initial_reporting_parameters$sd
    )
  }

  cl <- makeCluster(4)
  registerDoParallel(cl)
  res <- incidenceinflation::mcmc(niterations=niterations,
                                  df,
                                  priors,
                                  serial_parameters,
                                  initial_cases_true,
                                  initial_reporting_parameters,
                                  initial_Rt,
                                  reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                                  serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
                                  nchains=4, is_parallel=TRUE,
                                  is_negative_binomial=is_negative_binomial,
                                  initial_overdispersion=10,
                                  overdispersion_metropolis_sd=1)
  stopCluster(cl)
  
  res
}

fit_stan_model <- function(model,
                      df,
                      w_max,
                      generation_ws,
                      window,
                      I_init,
                      n_iterations,
                      include_opt_init=FALSE,
                      include_optimisation=TRUE) {
  options(mc.cores=4)
  
  data_stan <- list(
    N=nrow(df),
    K=max(window),
    window=window,
    C=df$C,
    wmax=w_max,
    w=generation_ws,
    I_history=I_init
  )
  
  
  if(!include_optimisation) {
    fit <- sampling(model, data=data_stan,
                    iter=n_iterations, chains=4)
    return(fit)
  }
  
  # optimise then sample
  if(include_opt_init) {
    
    init_fn <- function() {
      list(R=rep(1, max(window)))
    }
    fit <- optimizing(model, data=data_stan, as_vector=FALSE, init=init_fn)
    
  } else {
    fit <- optimizing(model, data=data_stan, as_vector=FALSE)
  }
  init_fn <- function() {
    list(R=fit$par$R,
         sigma=fit$par$sigma)
  }
  
  fit <- sampling(model, data=data_stan,
                  iter=n_iterations, chains=4,
                  init=init_fn)
  fit
}

fit_both_stan <- function(model,
                          df_full,
                          df_snapshot,
                          w_max,
                          generation_ws,
                          window,
                          I_init,
                          n_iterations,
                          include_opt_init=FALSE,
                          include_optimisation=TRUE) {
  
  fit_full <- fit_stan_model(model,
                             df_full,
                             w_max,
                             generation_ws,
                             window,
                             I_init,
                             n_iterations,
                             include_opt_init,
                             include_optimisation)
  fit_snapshot <- fit_stan_model(model,
                             df_snapshot,
                             w_max,
                             generation_ws,
                             window,
                             I_init,
                             n_iterations,
                             include_opt_init,
                             include_optimisation)
  list(full=fit_full,
       snapshot=fit_snapshot)
}

fit_model_uncertain_cases_stan <- function(df,
                      serial_parameters,
                      niterations=400,
                      nchains=1,
                      is_negative_binomial=FALSE,
                      is_maximise=FALSE,
                      initial_cases_true=NULL,
                      initial_reporting_parameters=NULL,
                      initial_Rt=NULL,
                      initial_overdispersion=10,
                      use_stan_sampling=FALSE,
                      n_stan_warmup=6,
                      n_stan_iterations=3,
                      step_size=NULL
                      ) {
  
  # number of reporting distributions
  num_reporting_pieces <- 1
  if("reporting_piece_index" %in% colnames(df)) {
    num_reporting_pieces <- max(df$reporting_piece_index)
  }
  
  # initial values of parameter estimates
  if(is.null(initial_cases_true)) {
    initial_cases_true=df %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=last(cases_reported))
  }
  
  if(is.null(initial_reporting_parameters)) {
    initial_reporting_parameters = list(location=10, scale=3)
  }
  
  if(is.null(initial_Rt)) {
    Rt_index_count <- max(df$Rt_index)
    initial_Rt <- tibble(
      Rt_index=seq(1, Rt_index_count, 1),
      Rt=1
    )
  }
  
  # priors
  Rt_prior <- list(location=5, scale=5, is_gamma=TRUE) # equates to gamma(shape=1, rate=0.2)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=10,
                                mean_sigma=15,
                                sd_mu=5,
                                sd_sigma=15), ## uninformative reporting prior
                 max_cases=max(df$cases_reported) * 10) ## wide range of cases
  
  if(is_negative_binomial) {
    print("Fitting using a NB renewal model...")
    priors$overdispersion <- list(location=0, scale=10)
  }
  
  if(num_reporting_pieces > 1) {
    initial_reporting_parameters <- tibble(
      reporting_piece_index=seq(1, num_reporting_pieces, 1),
      location=initial_reporting_parameters$location,
      scale=initial_reporting_parameters$scale
    )
  }
  
  stan_model <- cmdstanr::cmdstan_model("src/stan/conditional_renewal.stan",
                                        compile_model_methods=TRUE,
                                        force_recompile=TRUE)
  res <- incidenceinflationstan::mcmc(niterations=niterations,
              df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              stan_model,
              maximise=is_maximise,
              is_negative_binomial=is_negative_binomial,
              initial_overdispersion=initial_overdispersion,
              use_stan_sampling=use_stan_sampling,
              nchains=nchains,
              n_stan_warmup=n_stan_warmup,
              n_stan_iterations=n_stan_iterations,
              step_size=step_size)
  
  res
}