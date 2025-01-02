# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble", "dplyr", "ggplot2", "purrr", "forcats",
               "incidenceinflation", "incidenceinflationstan",
               "latex2exp",
               "doParallel", "foreach",
               "posterior", "rstan",
               "cowplot", "tidyr", "lubridate")
)

# options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
source("src/r/utils.R")
source("src/r/generate_data.R")
source("src/r/fit_models.R")
source("src/r/plotting.R")
source("src/r/mcmc_summary.R")
source("src/r/outbreak_over.R")
source("src/r/beta_binomial_noise.R")

# global parameters
n_iterations <- 200
niterations_step_down <- 3000
niterations_step_down_stan <- 2000
niterations_measles <- 4000
niterations_dengue <- 2000
niterations_resurgence <- 2000

# Replace the target list below with your own:
list(
  
  # targets used in publication
  tar_target(ultimate_targets, {
    list(
        file_plot_three_panel_step_down,
        file_plot_step_down_r_cases,
        file_plot_cases_dengue_combined,
        file_plot_dengue_r_estimates,
        file_plot_cases_measles_combined,
        file_plot_measles_r_estimates
         # file_plot_ebola_real_rt_noisy,
         # file_plot_ebola_real_rt_noisy_all,
         # file_plot_ebola_real_rt,
         # file_plot_ebola_real_rt_all,
         # file_plot_cases_ebola_noisy_combined,
         # file_plot_cases_ebola_combined
     )
  }),
  
  # generate simulated covid dataset
  tar_target(serial_parameters_covid, {
    # Picked as gamma distribution close to best fit
    # distribution from Nishiura et al. 2020
    mu <- 4.6
    sigma <- 4.8
    list(mean=mu, sd=sigma)
  }),
  tar_target(rt_step_down, {
    v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
    v_Rt
  }),
  tar_target(rt_step_down_df, {
    tibble(Rt=rt_step_down,
           time_onset=seq_along(Rt))
  }),
  tar_target(reporting_parameters_10_3, {
    list(mean=10, sd=3)
  }),
  tar_target(data_sim_step_down_up,
             generate_data(rt_step_down,
                           serial_parameters_covid,
                           reporting_parameters_10_3)),
  tar_target(serial_parameters_ebola, {
    # taken from all countries 2015 estimate from Van Kerkhove et al.
    # 2014 (their table 6)
    mu <- 15.3
    sigma <- 9.3
    
    list(mean=mu, sd=sigma)}
  ),
  
  # fit ebola dataset across a range of "current" times
  tar_target(results_mcmc_resurgence_full,
             fit_ebola(0, data_sim_resurgence_ebola,
                       serial_parameters_ebola,
                       niterations_resurgence)),
  tar_target(results_mcmc_resurgence_25,
             fit_ebola(25, data_sim_resurgence_ebola,
                       serial_parameters_ebola,
                       niterations_resurgence)),
  tar_target(results_mcmc_resurgence_50,
             fit_ebola(50, data_sim_resurgence_ebola,
                       serial_parameters_ebola,
                       niterations_resurgence)),
  tar_target(results_mcmc_resurgence_75,
             fit_ebola(75, data_sim_resurgence_ebola,
                       serial_parameters_ebola,
                       niterations_resurgence)),
  tar_target(summary_mcmc_resurgence_full,
             mcmc_summary(results_mcmc_resurgence_full)),
  tar_target(summary_mcmc_resurgence_25,
             mcmc_summary(results_mcmc_resurgence_25)),
  tar_target(summary_mcmc_resurgence_50,
             mcmc_summary(results_mcmc_resurgence_50)),
  tar_target(summary_mcmc_resurgence_75,
             mcmc_summary(results_mcmc_resurgence_75)),
  tar_target(cases_resurgence_combined, {
    
    df_full <- results_mcmc_resurgence_full$cases %>% 
      mutate(max_time="125")
    df_100 <- results_mcmc_resurgence_25$cases %>% 
      mutate(max_time="100")
    df_75 <- results_mcmc_resurgence_50$cases %>% 
      mutate(max_time="75")
    df_50 <- results_mcmc_resurgence_75$cases %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_100,
        df_75,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "75", "100", "125"))
  }),
  tar_target(plot_cases_resurgence_combined,
             plot_cases_vs_true_series_1(cases_resurgence_combined,
                                         data_sim_resurgence_ebola,
                                         include_reporting=TRUE)),
  tar_target(file_plot_cases_resurgence_combined, {
    filename <- "figures/resurgence_cases.pdf"
    ggsave(filename,
           plot_cases_resurgence_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(rt_resurgence_combined, {
    
    df_full <- results_mcmc_resurgence_full$Rt %>% 
      mutate(max_time="125")
    df_100 <- results_mcmc_resurgence_25$Rt %>% 
      mutate(max_time="100")
    df_75 <- results_mcmc_resurgence_50$Rt %>% 
      mutate(max_time="75")
    df_50 <- results_mcmc_resurgence_75$Rt %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_100,
        df_75,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "75", "100", "125"))
  }),
  tar_target(plot_rt_resurgence_combined,
             plot_rt_vs_true_series(rt_resurgence_combined,
                                    resurgence_rt,
                                    data_sim_resurgence_ebola,
                                    threshold_onset_time=21)),
  tar_target(file_plot_rt_resurgence_combined, {
    filename <- "figures/resurgence_rt.pdf"
    ggsave(filename,
           plot_rt_resurgence_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(reporting_delay_resurgence_combined, {
    
    df_full <- results_mcmc_resurgence_full$reporting %>% 
      mutate(max_time="125")
    df_100 <- results_mcmc_resurgence_25$reporting %>% 
      mutate(max_time="100")
    df_75 <- results_mcmc_resurgence_50$reporting %>% 
      mutate(max_time="75")
    df_50 <- results_mcmc_resurgence_75$reporting %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_100,
        df_75,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "75", "100", "125"))
  }),
  tar_target(plot_reporting_delay_resurgence_combined,
             plot_delay_true_series(reporting_delay_resurgence_combined,
                                    reporting_parameters_8_8)),
  tar_target(file_plot_reporting_delay_resurgence_combined, {
    filename <- "figures/resurgence_reporting_delays.pdf"
    ggsave(filename,
           plot_reporting_delay_resurgence_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(model_nb_raw, "src/stan/negative_binomial.stan",
             format = "file"),
  tar_target(model_nb, {
    rstan_options(auto_write=TRUE)
    stan_model(model_nb_raw)
  }),
  tar_target(ebola_processed_simple, {
    data_sim_resurgence_ebola %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported))
  }),
  tar_target(ebola_generation_ws,
             weights_series(40, serial_parameters_ebola)),
  tar_target(ebola_rt_index, {
    data_sim_resurgence_ebola %>% 
      select(time_onset, Rt_index) %>% 
      unique()
  }),
  tar_target(fits_ebola_stan_0, fit_both_stan_simple(
    ebola_processed_simple, data_sim_resurgence_ebola,
    ebola_generation_ws,
    ebola_rt_index, nrow(ebola_rt_index), niterations_resurgence, model_nb)),
  tar_target(fits_ebola_stan_25, fit_both_stan_simple(
    ebola_processed_simple, data_sim_resurgence_ebola,
    ebola_generation_ws,
    ebola_rt_index, nrow(ebola_rt_index) - 25, niterations_resurgence, model_nb)),
  tar_target(fits_ebola_stan_50, fit_both_stan_simple(
    ebola_processed_simple, data_sim_resurgence_ebola,
    ebola_generation_ws,
    ebola_rt_index, nrow(ebola_rt_index) - 50, niterations_resurgence, model_nb)),
  tar_target(fits_ebola_stan_75, fit_both_stan_simple(
    ebola_processed_simple, data_sim_resurgence_ebola,
    ebola_generation_ws,
    ebola_rt_index, nrow(ebola_rt_index) - 75, niterations_resurgence, model_nb)),
  tar_target(summarise_ebola_r_estimates, {
    
    # process stan fits
    R_0 <- summarise_r(fits_ebola_stan_0) %>% 
      mutate(max_time=0)
    R_25 <- summarise_r(fits_ebola_stan_25) %>% 
      mutate(max_time=25)
    R_50 <- summarise_r(fits_ebola_stan_50) %>% 
      mutate(max_time=50)
    R_75 <- summarise_r(fits_ebola_stan_75) %>% 
      mutate(max_time=75)
    R_stan <- R_0 %>% 
      bind_rows(
        R_0,
        R_25,
        R_50,
        R_75
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_0 <- quantiles_r_other(results_mcmc_resurgence_full) %>% 
      mutate(max_time=0)
    R_25 <- quantiles_r_other(results_mcmc_resurgence_25) %>% 
      mutate(max_time=25)
    R_50 <- quantiles_r_other(results_mcmc_resurgence_50) %>% 
      mutate(max_time=50)
    R_75 <- quantiles_r_other(results_mcmc_resurgence_75) %>% 
      mutate(max_time=75)
    R_inflation <- R_0 %>% 
      bind_rows(
        R_25,
        R_50,
        R_75
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(plot_ebola_r_estimates, {
    ebola_rt_index %>% 
      left_join(summarise_ebola_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      mutate(max_time=fct_relevel(max_time, "75", "50", "25", "0")) %>% 
      filter(time_onset >= 21) %>% 
      ggplot(aes(x=time_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_brewer("Estimator", palette = "Dark2") +
      scale_fill_brewer("Estimator",palette = "Dark2") +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~max_time)
  }),
  tar_target(file_plot_ebola_r_estimates, {
    filename <- "figures/resurgence_rt_comparison.pdf"
    ggsave(filename,
           plot_ebola_r_estimates,
           width = 10, height = 4);
    filename
  }),
  
  
  # fit covid dataset across a range of "current" times
  tar_target(results_mcmc_step_down_full,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  tar_target(results_mcmc_step_down_90,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true) %>% 
                         filter(time_reported <=90),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  tar_target(results_mcmc_step_down_80,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true) %>% 
                         filter(time_reported <=80),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  tar_target(results_mcmc_step_down_70,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true) %>% 
                         filter(time_reported <=70),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  tar_target(results_mcmc_step_down_60,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true) %>% 
                         filter(time_reported <=60),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  tar_target(results_mcmc_step_down_50,
             fit_model(data_sim_step_down_up %>% 
                         select(-cases_true) %>% 
                         filter(time_reported <=50),
                       serial_parameters_covid,
                       niterations=niterations_step_down)),
  
  
  # check convergence diagnostics
  tar_target(summary_mcmc_step_down_full,
             mcmc_summary(results_mcmc_step_down_full)),
  tar_target(summary_mcmc_step_down_90,
             mcmc_summary(results_mcmc_step_down_90)),
  tar_target(summary_mcmc_step_down_80,
             mcmc_summary(results_mcmc_step_down_80)),
  tar_target(summary_mcmc_step_down_70,
             mcmc_summary(results_mcmc_step_down_70)),
  tar_target(summary_mcmc_step_down_60,
             mcmc_summary(results_mcmc_step_down_60)),
  tar_target(summary_mcmc_step_down_50,
             mcmc_summary(results_mcmc_step_down_50)),
  tar_target(summary_mcmc_step_down_all, {
    list(
      summary_mcmc_step_down_full,
      summary_mcmc_step_down_90,
      summary_mcmc_step_down_80,
      summary_mcmc_step_down_70,
      summary_mcmc_step_down_60,
      summary_mcmc_step_down_50
    )
  }),
  
  ## plot covid results
  tar_target(cases_step_down_combined, {
    
    df_full <- results_mcmc_step_down_full$cases %>% 
      mutate(max_time="100")
    df_90 <- results_mcmc_step_down_90$cases %>% 
      mutate(max_time="90")
    df_80 <- results_mcmc_step_down_80$cases %>% 
      mutate(max_time="80")
    df_70 <- results_mcmc_step_down_70$cases %>% 
      mutate(max_time="70")
    df_60 <- results_mcmc_step_down_60$cases %>% 
      mutate(max_time="60")
    df_50 <- results_mcmc_step_down_50$cases %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_90,
        df_80,
        df_70,
        df_60,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "60", "70", "80", "90", "100"))
  }),
  tar_target(plot_cases_step_down_combined,
             plot_cases_vs_true_series_1(cases_step_down_combined,
                                       data_sim_step_down_up,
                                       include_reporting=TRUE)),
  tar_target(file_plot_cases_step_down_combined, {
    filename <- "figures/step_down_cases.pdf"
    ggsave(filename,
           plot_cases_step_down_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(rt_step_down_combined, {
    
    df_full <- results_mcmc_step_down_full$Rt %>% 
      mutate(max_time="100")
    df_90 <- results_mcmc_step_down_90$Rt %>% 
      mutate(max_time="90")
    df_80 <- results_mcmc_step_down_80$Rt %>% 
      mutate(max_time="80")
    df_70 <- results_mcmc_step_down_70$Rt %>% 
      mutate(max_time="70")
    df_60 <- results_mcmc_step_down_60$Rt %>% 
      mutate(max_time="60")
    df_50 <- results_mcmc_step_down_50$Rt %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_90,
        df_80,
        df_70,
        df_60,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "60", "70", "80", "90", "100"))
  }),
  tar_target(plot_rt_step_down_combined,
             plot_rt_vs_true_series(rt_step_down_combined,
                                    rt_step_down,
                                    data_sim_step_down_up,
                                    threshold_onset_time=21)),
  tar_target(file_plot_rt_step_down_combined, {
    filename <- "figures/step_down_rt.pdf"
    ggsave(filename,
           plot_rt_step_down_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(reporting_delay_step_down_combined, {
    
    df_full <- results_mcmc_step_down_full$reporting %>% 
      mutate(max_time="100")
    df_90 <- results_mcmc_step_down_90$reporting %>% 
      mutate(max_time="90")
    df_80 <- results_mcmc_step_down_80$reporting %>% 
      mutate(max_time="80")
    df_70 <- results_mcmc_step_down_70$reporting %>% 
      mutate(max_time="70")
    df_60 <- results_mcmc_step_down_60$reporting %>% 
      mutate(max_time="60")
    df_50 <- results_mcmc_step_down_50$reporting %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_90,
        df_80,
        df_70,
        df_60,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "60", "70", "80", "90", "100"))
  }),
  tar_target(plot_reporting_delay_step_down_combined,
             plot_delay_true_series(reporting_delay_step_down_combined,
                                    reporting_parameters_10_3)),
  tar_target(file_plot_reporting_delay_step_down_combined, {
    filename <- "figures/step_down_reporting_delays.pdf"
    ggsave(filename,
           plot_reporting_delay_step_down_combined,
           width = 10, height = 4);
    filename
  }),
  
  ## three panel figure
  tar_target(plot_cases_step_down_simple,
             plot_cases_simple(results_mcmc_step_down_full$cases %>% 
                                 mutate(max_time="100"),
                               data_sim_step_down_up %>% filter(time_onset>=20),
                               include_reporting=TRUE)),
  tar_target(plot_reporting_delay_step_down_simple,
             plot_delay_true_series_simple(results_mcmc_step_down_full$reporting %>% 
                                             mutate(max_time="100"),
                                           reporting_parameters_10_3)),
  tar_target(step_down_processed_simple, {
    data_sim_step_down_up %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_true))
  }),
  tar_target(step_down_processed_snapshot, {
    data_sim_step_down_up %>% 
      filter(time_onset <= max(time_onset)) %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported))
  }),
  tar_target(step_down_generation_ws,
             weights_series(40, serial_parameters_covid)),
  tar_target(step_down_rt_index, {
    data_sim_step_down_up %>% 
      select(time_onset, Rt_index) %>% 
      unique()
  }),
  tar_target(fits_step_down_stan_full, {
    fit_both_stan(
      model_poisson,
      step_down_processed_simple,
      step_down_processed_snapshot,
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index,
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(fits_step_down_stan_90, {
    
    final_time <- 100 - 10
    
    fit_both_stan(
      model_poisson,
      data_sim_step_down_up %>% 
        filter(time_onset <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_true)),
      data_sim_step_down_up %>% 
        filter(time_reported <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index[1:final_time],
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(fits_step_down_stan_80, {
    
    final_time <- 100 - 20
    
    fit_both_stan(
      model_poisson,
      data_sim_step_down_up %>% 
        filter(time_onset <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_true)),
      data_sim_step_down_up %>% 
        filter(time_reported <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index[1:final_time],
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(fits_step_down_stan_70, {
    
    final_time <- 100 - 30
    
    fit_both_stan(
      model_poisson,
      data_sim_step_down_up %>% 
        filter(time_onset <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_true)),
      data_sim_step_down_up %>% 
        filter(time_reported <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index[1:final_time],
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(fits_step_down_stan_60, {
    
    final_time <- 100 - 40
    
    fit_both_stan(
      model_poisson,
      data_sim_step_down_up %>% 
        filter(time_onset <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_true)),
      data_sim_step_down_up %>% 
        filter(time_reported <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index[1:final_time],
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(fits_step_down_stan_50, {
    
    final_time <- 100 - 50
    
    fit_both_stan(
      model_poisson,
      data_sim_step_down_up %>% 
        filter(time_onset <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_true)),
      data_sim_step_down_up %>% 
        filter(time_reported <= final_time) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      step_down_generation_ws,
      step_down_rt_index$Rt_index[1:final_time],
      I_init=rep(1, 40),
      n_iterations=niterations_step_down_stan)
  }),
  tar_target(summarise_step_down_r_estimates, {
    
    # process stan fits
    R_0 <- summarise_r(fits_step_down_stan_full) %>% 
      mutate(max_time=100)
    R_10 <- summarise_r(fits_step_down_stan_90) %>% 
      mutate(max_time=90)
    R_20 <- summarise_r(fits_step_down_stan_80) %>% 
      mutate(max_time=80)
    R_30 <- summarise_r(fits_step_down_stan_70) %>% 
      mutate(max_time=70)
    R_40 <- summarise_r(fits_step_down_stan_60) %>% 
      mutate(max_time=60)
    R_50 <- summarise_r(fits_step_down_stan_50) %>% 
      mutate(max_time=50)
    R_stan <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_0 <- quantiles_r_other(results_mcmc_step_down_full) %>% 
      mutate(max_time=100)
    R_10 <- quantiles_r_other(results_mcmc_step_down_90) %>% 
      mutate(max_time=90)
    R_20 <- quantiles_r_other(results_mcmc_step_down_80) %>% 
      mutate(max_time=80)
    R_30 <- quantiles_r_other(results_mcmc_step_down_70) %>% 
      mutate(max_time=70)
    R_40 <- quantiles_r_other(results_mcmc_step_down_60) %>% 
      mutate(max_time=60)
    R_50 <- quantiles_r_other(results_mcmc_step_down_50) %>% 
      mutate(max_time=50)
    R_inflation <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50
      ) %>% 
      mutate(estimator="incidenceinflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(plot_step_down_r_estimates_full, {
    
    df <- step_down_rt_index %>% 
      left_join(summarise_step_down_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      mutate(variant=paste0(estimator, "-", type)) %>% 
      filter(max_time==100) %>% 
      filter(time_onset >= 20) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="incidenceinflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
      df %>% 
      filter(time_onset>20) %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper, y=middle, fill=variant), alpha=0.4) +
      geom_line(aes(colour=variant, y=middle, fill=variant)) +
      geom_line(data=rt_step_down_df,
                aes(y=Rt), linetype=2) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      xlim(20, NA) +
      theme(legend.position = c(0.65, 0.4))
  }),
  tar_target(plot_three_panel_step_down, {
    plot_grid(plot_cases_step_down_simple + xlim(20, NA),
              plot_step_down_r_estimates_full,
              plot_reporting_delay_step_down_simple,
              labels = c("A.", "B.", "C."),
              nrow=3)
  }),
  tar_target(file_plot_three_panel_step_down, {
    filename <- "figures/step_down_three_panel.pdf"
    ggsave(filename,
           plot_three_panel_step_down,
           width = 10, height = 8);
    filename
  }),
  tar_target(file_plot_three_panel_step_down_png, {
    filename <- "figures/step_down_three_panel.png"
    ggsave(filename,
           plot_three_panel_step_down,
           width = 10, height = 8, dpi = 1000);
    filename
  }),
  
  tar_target(plot_step_down_r_estimates, {
    df <- step_down_rt_index %>% 
      left_join(summarise_step_down_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "60", "70", "80", "90", "100")) %>% 
      filter(time_onset >= 20) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="incidenceinflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    df %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=variant), alpha=0.4) +
      geom_line(aes(colour=variant, y=middle)) +
      geom_line(data=rt_step_down_df %>% filter(time_onset >= 20),
                aes(y=Rt), linetype=2) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset time") +
      facet_wrap(~max_time) +
      theme(
        legend.position = "bottom"
      )
  }),
  tar_target(file_plot_step_down_rt_estimates, {
    filename <- "figures/step_down_rt_all_estimates.pdf"
    ggsave(filename,
           plot_step_down_r_estimates,
           width = 10, height = 4);
    filename
  }),
  tar_target(plot_step_down_r_cases, {
    plot_grid(plot_cases_step_down_combined,
              plot_step_down_r_estimates,
              labels = c("A.", "B."),
              nrow=2)
  }),
  tar_target(file_plot_step_down_r_cases, {
    filename <- "figures/step_down_rt_rt_cases.pdf"
    ggsave(filename,
           plot_step_down_r_cases,
           width = 10, height = 8);
    filename
  }),
  tar_target(file_plot_step_down_r_cases_png, {
    filename <- "figures/step_down_rt_rt_cases.png"
    ggsave(filename,
           plot_step_down_r_cases,
           width = 10, height = 8, dpi=1000);
    filename
  }),
  
  # ebola-like outbreak with shortening delays
  tar_target(reporting_parameters_8_8, {
    tribble(mean=8, sd=8)
  }),
  tar_target(resurgence_rt_short, {
    R1 <- 2.5
    R2 <- 0.5
    R3 <- 2.5
    v_R <- c(rep(R1, 30), rep(R2, 30), rep(R3, 30))
    v_R
  }),
  tar_target(resurgence_rt_df, {
    tibble(Rt=resurgence_rt_short,
           time_onset=seq_along(Rt))
  }),
  tar_target(reporting_parameters_shortening, {
    indices <-tribble(
      ~reporting_piece_index, ~mean, ~sd,
      1, 12, 5,
      2, 5, 3)
    reporting_piece_indices <- unlist(map(seq(1, 2, 1), ~rep(., 45)))
    lookup <- tibble(time_onset=seq_along(resurgence_rt_short),
                     reporting_piece_index=reporting_piece_indices)
    indices %>% 
      left_join(lookup, by="reporting_piece_index")
  }),
  tar_target(data_sim_shortening_ebola,
             generate_data(resurgence_rt_short,
                           serial_parameters_ebola,
                           reporting_parameters_shortening,
                           kappa=5,
                           days_total=length(resurgence_rt_short))),
  tar_target(shortening_ebola_rt_index, {
    data_sim_shortening_ebola %>% 
      select(time_onset, Rt_index) %>% 
      unique()
  }),
  tar_target(opt_shortening_ebola_full, {
    fit_model_uncertain_cases_stan(
      data_sim_shortening_ebola %>% 
      thin_series() %>% 
      left_join(reporting_parameters_shortening, by="time_onset") %>% 
      select(-c(mean, sd, cases_true)),
      serial_parameters_ebola,
      niterations=10,
      is_negative_binomial=TRUE,
      is_maximise=TRUE)
  }),
  tar_target(fit_shortening_ebola_full, {
    
    initial_cases_true <- opt_shortening_ebola_full$cases %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_Rt <- opt_shortening_ebola_full$Rt %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_reporting_parameters <- opt_shortening_ebola_full$reporting %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_overdispersion <- opt_shortening_ebola_full$other %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration) %>% 
      pull(overdispersion)
    
    fit_model_uncertain_cases_stan(
      data_sim_shortening_ebola %>% 
        thin_series() %>% 
        left_join(reporting_parameters_shortening, by="time_onset") %>% 
        select(-c(mean, sd, cases_true)),
      serial_parameters_ebola,
      niterations=100,
      is_negative_binomial=TRUE,
      is_maximise=FALSE,
      initial_cases_true=initial_cases_true,
      initial_reporting_parameters=initial_reporting_parameters,
      initial_Rt=initial_Rt,
      initial_overdispersion=initial_overdispersion,
      use_stan_sampling=TRUE,
      step_size=0.5e-3,
      nchains=2)
  }),
  tar_target(summary_fit_shortening_ebola_full,
             mcmc_summary_stan_uncertain_cases(fit_shortening_ebola_full)),
  tar_target(summarise_shortening_ebola_full_r_estimates, {
    
    quantiles_r_other(fit_shortening_ebola_full) %>% 
      mutate(estimator="inflation-snapshot")
  }),
  tar_target(plot_rt_shortening_ebola_full,
             {
               shortening_ebola_rt_index %>% 
                 left_join(summarise_shortening_ebola_full_r_estimates,
                           by="Rt_index",
                           relationship =
                             "many-to-many") %>%
                 left_join(resurgence_rt_df, by="time_onset") %>% 
                 mutate(variant=paste0(estimator, "-", type)) %>% 
                 filter(time_onset >= 20) %>% 
                 ggplot(aes(x=time_onset, y=middle)) +
                 geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4,
                             fill="#1B9E77") +
                 geom_line(colour="#1B9E77") +
                 geom_line(aes(y=Rt), linetype=2) +
                 ylab(TeX("$R_t$")) +
                 xlab("Onset time") +
                 theme(legend.position="none") +
                 xlim(20, 90)
             }),
  tar_target(plot_shortening_ebola_full_delay, {
    data_sim_shortening_ebola %>% 
      mutate(time_delay=time_reported-time_onset) %>%
      ggplot(aes(x=time_onset, y=time_delay,
                 fill=cases_reported)) +
      geom_tile() +
      scale_fill_distiller("Cases reported", direction=1) +
      xlab("Onset date") +
      ylab("Reporting delay") +
      theme(
        text=element_text(size=16)
      )
  }),
  tar_target(plot_shortening_ebola_full_delay_1, {
    
    cols <- c("long"="#E7298A", "short"="#66A61E")
    data_sim_shortening_ebola %>%
      mutate(is_short_delay=if_else(time_onset > 50, "short", "long")) %>%
      filter(time_onset >= 20) %>% 
      ggplot(aes(x=time_reported,
                 y=cases_reported,
                 group=as.factor(time_onset),
                 colour=is_short_delay)) +
      geom_line() +
      scale_color_manual("Rep. delay", values = cols) +
      xlab("Onset time & report time") +
      ylab("Cases reported") +
      scale_y_sqrt(limits=c(0, 3000)) +
      theme(legend.position = c(0.7, 0.75))
  }),
  tar_target(plot_shortening_ebola_full_cases, {
    
    df_observed <- data_sim_shortening_ebola %>% 
      group_by(time_onset) %>% 
      summarise(
        cases_observed=last(cases_reported)
      )
    
    fit_shortening_ebola_full$cases %>% 
      group_by(time_onset) %>% 
      summarise(
        middle=median(cases_true),
        lower=quantile(cases_true, 0.975),
        upper=quantile(cases_true, 0.025)
        ) %>% 
      left_join(data_sim_shortening_ebola, by="time_onset") %>% 
      left_join(df_observed, by="time_onset") %>% 
      filter(time_onset >= 20) %>% 
      ggplot(aes(x=time_onset)) +
      geom_ribbon(aes(ymin=lower, ymax=upper),
                  fill="#1B9E77", alpha=0.3) +
      geom_line(aes(y=middle), colour="#1B9E77") +
      geom_line(aes(y=cases_observed), colour="orange") +
      geom_line(aes(y=cases_true), linetype=2) +
      scale_y_sqrt() +
      xlab("Onset time & report time") +
      ylab("Cases")
  }),
  tar_target(plot_reporting_delay_shortening_ebola_full, {
    
    rep_df <- fit_shortening_ebola_full$reporting
    
    x <- seq(0, 25, 0.1)
    for(i in seq_along(rep_df$reporting_piece_index)) {
      
      mu <- rep_df$location[i]
      sigma <- rep_df$scale[i]
      a <- mu ^ 2 / (sigma ^ 2)
      b <- mu / (sigma ^ 2)
      
      dens <- dgamma(x, a, b)
      temp_df <- tibble(x=x, pdf=dens,
                        iteration=rep_df$iteration[i],
                        reporting_piece_index=rep_df$reporting_piece_index[i])
      if(i == 1)
        big_df <- temp_df
      else
        big_df <- big_df %>% bind_rows(temp_df)
    }
    
    true_mus <- c(12, 5)
    true_sds <- c(5, 3)
    for(i in seq_along(true_mus)) {
      
      mu <- true_mus[i]
      sigma <- true_sds[i]
      a <- mu ^ 2 / (sigma ^ 2)
      b <- mu / (sigma ^ 2)
      
      dens <- dgamma(x, a, b)
      temp_df <- tibble(x=x, pdf=dens,
                        reporting_piece_index=i)
      if(i == 1)
        true_df <- temp_df
      else
        true_df <- true_df %>% bind_rows(temp_df)
    }
    
    cols <- c("1"="#E7298A", "2"="#66A61E")
    big_df %>% 
      mutate(comb=paste0(iteration, reporting_piece_index)) %>% 
      mutate(comb_fac=as.factor(comb)) %>% 
      ggplot(aes(x=x, y=pdf)) +
      geom_line(alpha=0.5, aes(group=as.factor(comb_fac),
                               colour=as.factor(reporting_piece_index))) +
      geom_line(data=true_df, linetype=2,
                aes(group=as.factor(reporting_piece_index))) +
      scale_colour_manual("Rep. delay", values=cols) +
      xlab("Reporting delay, days") +
      ylab("Probability\ndensity") +
      theme(legend.position = "none",
            text = element_text(size=8))
  }),
  tar_target(plot_ebola_full_delay_combined, {
    plot_shortening_ebola_full_delay_1 +
      annotation_custom(ggplotGrob(plot_reporting_delay_shortening_ebola_full), 
                        xmin = 17, xmax = 35, ymin = 16, ymax = 56)
  }),
  tar_target(plot_three_panel_shortening_ebola, {
    plot_grid(plot_ebola_full_delay_combined,
              plot_shortening_ebola_full_cases,
              plot_rt_shortening_ebola_full,
              labels = c("A.", "B.", "C."),
              nrow=3)
  }),
  tar_target(file_plot_three_panel_shortening_ebola, {
    filename <- "figures/ebola_shortening_three_panel.pdf"
    ggsave(filename,
           plot_three_panel_shortening_ebola,
           width = 12, height = 8);
    filename
  }),
  tar_target(file_plot_three_panel_shortening_ebola_png, {
    filename <- "figures/ebola_shortening_three_panel.png"
    ggsave(filename,
           plot_three_panel_shortening_ebola,
           width = 12, height = 8, dpi=1000);
    filename
  }),
  
  ## try now allowing a reporting delay for different subperiods
  tar_target(reporting_parameters_shortening_one_value, {
    reporting_parameters_shortening %>% 
      mutate(reporting_piece_index=1)
  }),
  tar_target(opt_shortening_ebola_full_one_value, {
    fit_model_uncertain_cases_stan(
      data_sim_shortening_ebola %>% 
        thin_series() %>% 
        left_join(reporting_parameters_shortening_one_value, by="time_onset") %>% 
        select(-c(mean, sd, cases_true)),
      serial_parameters_ebola,
      niterations=15,
      is_negative_binomial=TRUE,
      is_maximise=TRUE)
  }),
  tar_target(reporting_parameters_shortening_lots_of_vals, {
    
    indices <- unlist(map(seq(1, 3, 1), ~rep(., 30)))
    
    reporting_parameters_shortening %>% 
      mutate(reporting_piece_index=indices)
  }),
  tar_target(opt_shortening_ebola_full_lots_of_vals, {
    fit_model_uncertain_cases_stan(
      data_sim_shortening_ebola %>% 
        thin_series() %>% 
        left_join(reporting_parameters_shortening_lots_of_vals, by="time_onset") %>% 
        select(-c(mean, sd, cases_true)),
      serial_parameters_ebola,
      niterations=10,
      is_negative_binomial=TRUE,
      is_maximise=TRUE)
  }),
  tar_target(reporting_delay_progress_lots_of_vals, {
    
    rep_df <- opt_shortening_ebola_full_lots_of_vals$reporting %>% 
      mutate(
        lower=location-scale,
        upper=location+scale
      ) %>% 
      mutate(period=case_when(
        reporting_piece_index==1~"<=30 days",
        reporting_piece_index==2~"31-60 days",
        reporting_piece_index==3~">=61 days"
      )) %>% 
      mutate(period=as.factor(period)) %>% 
      mutate(period=fct_relevel(period, "<=30 days", "31-60 days", ">=61 days"))
    true_df <- tribble(
      ~reporting_piece_index, ~location, ~lower, ~upper,
      1, 12, 7, 17,
      2, 1/2 * 12 + 1/2 * 5, NA, NA,
      3, 5, 2, 8
    )
    short_df <- rep_df %>% 
      select(reporting_piece_index, iteration) %>% 
      unique() %>% 
      left_join(true_df, by="reporting_piece_index")%>% 
      mutate(period=case_when(
        reporting_piece_index==1~"<=30 days",
        reporting_piece_index==2~"31-60 days",
        reporting_piece_index==3~">=61 days"
      )) %>% 
      mutate(period=as.factor(period)) %>% 
      mutate(period=fct_relevel(period, "<=30 days", "31-60 days", ">=61 days"))
     rep_df %>%
       ggplot(aes(x=iteration, y=location)) +
       geom_line(colour="blue") +
       geom_ribbon(aes(ymin=lower, ymax=upper),
                   fill="blue", alpha=0.4) +
       geom_line(data=short_df, linetype=2) +
       geom_ribbon(data=short_df, aes(ymin=lower, ymax=upper),
                   fill="black", alpha=0.4) +
       scale_x_reverse(breaks=seq(1, 10, 1), position = "top") +
       scale_y_continuous(position = "right") +
       coord_flip() +
       facet_wrap(~period) +
       xlab("Iteration") +
       ylab("Reporting delay, days")
  }),
  tar_target(reporting_delay_progress_one_value, {
    
    rep_df <- opt_shortening_ebola_full_one_value$reporting %>% 
      mutate(
        lower=location-scale,
        upper=location+scale
      ) %>% 
      mutate(period="1-90 days") %>% 
      filter(iteration <= 10)
    true_df <- tribble(
      ~reporting_piece_index, ~location,
      1, 12 * 5/9 + 5 * 4/9
    )
    short_df <- rep_df %>% 
      select(reporting_piece_index, iteration) %>% 
      unique() %>% 
      left_join(true_df, by="reporting_piece_index")%>% 
      mutate(period=case_when(
        reporting_piece_index==1~"1-90 days",
      ))
    rep_df %>%
      ggplot(aes(x=iteration, y=location)) +
      geom_line(colour="blue") +
      geom_ribbon(aes(ymin=lower, ymax=upper),
                  fill="blue", alpha=0.4) +
      geom_line(data=short_df, linetype=2) +
      scale_x_reverse(breaks=seq(1, 10, 1), position = "top") +
      scale_y_continuous(position = "right") +
      coord_flip() +
      facet_wrap(~period) +
      xlab("Iteration") +
      ylab("Reporting delay, days")
  }),
  tar_target(reporting_delay_progress_two_vals, {
    
    rep_df <- opt_shortening_ebola_full$reporting %>% 
      mutate(
        lower=location-scale,
        upper=location+scale
      ) %>% 
      mutate(period=case_when(
        reporting_piece_index==1~"<=50 days",
        reporting_piece_index==2~">50 days"
      )) %>% 
      mutate(period=as.factor(period)) %>% 
      mutate(period=fct_relevel(period, "<=50 days", ">50 days"))
    true_df <- tribble(
      ~reporting_piece_index, ~location, ~lower, ~upper,
      1, 12, 7, 17,
      2, 5, 2, 8
    )
    short_df <- rep_df %>% 
      select(reporting_piece_index, iteration) %>% 
      unique() %>% 
      left_join(true_df, by="reporting_piece_index")%>% 
      mutate(period=case_when(
        reporting_piece_index==1~"<=50 days",
        reporting_piece_index==2~">50 days"
      )) %>% 
      mutate(period=as.factor(period)) %>% 
      mutate(period=fct_relevel(period, "<=50 days", ">50 days"))
    rep_df %>%
      ggplot(aes(x=iteration, y=location)) +
      geom_line(colour="blue") +
      geom_ribbon(aes(ymin=lower, ymax=upper),
                  fill="blue", alpha=0.4) +
      geom_ribbon(data=short_df, aes(ymin=lower, ymax=upper),
                  fill="black", alpha=0.4) +
      geom_line(data=short_df, linetype=2) +
      scale_x_reverse(breaks=seq(1, 10, 1), position = "top") +
      scale_y_continuous(position = "right") +
      coord_flip() +
      facet_wrap(~period) +
      xlab("Iteration") +
      ylab("Reporting delay, days")
  }),
  tar_target(reporting_delay_progress_stacked, {
    plot_grid(reporting_delay_progress_lots_of_vals,
              reporting_delay_progress_two_vals,
              reporting_delay_progress_one_value,
              nrow=3,
              labels = c("A.", "B.", "C."))
  }),
  tar_target(file_reporting_delay_progress_stacked, {
    
    filename <- "figures/ebola_reporting_delays_progress.pdf"
    ggsave(filename,
           reporting_delay_progress_stacked,
           width = 10, height = 8);
    filename
  }),
  tar_target(file_reporting_delay_progress_stacked_png, {
    
    filename <- "figures/ebola_reporting_delays_progress.png"
    ggsave(filename,
           reporting_delay_progress_stacked,
           width = 10, height = 8, dpi=1000);
    filename
  }),
  tar_target(plot_cases_reporting_delay_progress, {
    
    start_point <- 75
    cases_df_2 <- opt_shortening_ebola_full$cases %>% 
      mutate(period="B. days 1-50 & 51-90")
    cases_df_1 <- opt_shortening_ebola_full_one_value$cases %>% 
      mutate(period="A. days 1-90") %>% 
      filter(iteration <= 10)
    cases_df_3 <- opt_shortening_ebola_full_lots_of_vals$cases %>% 
      mutate(period="C. days 1-30 & 31-60 & 61-90")
    cases_df <- cases_df_1 %>% 
      bind_rows(
        cases_df_2,
        cases_df_3
      ) %>% 
      mutate(period=as.factor(period)) %>% 
      mutate(period=fct_relevel(period, "A. days 1-90", "B. days 1-50 & 51-90",
                                "C. days 1-30 & 31-60 & 61-90"))
    true_df <- data_sim_shortening_ebola %>% 
      group_by(time_onset) %>% 
      summarise(
        cases_true=first(cases_true)
      ) %>% 
      filter(time_onset >= start_point)
    observed_df <- data_sim_shortening_ebola %>% 
      group_by(time_onset) %>% 
      summarise(
        cases_true=last(cases_reported)
      ) %>% 
      filter(time_onset >= start_point)
    
    cases_df %>% 
      filter(time_onset >= start_point) %>% 
      ggplot(aes(x=time_onset, y=cases_true)) +
      geom_line(aes(colour=iteration,
                    group=as.factor(iteration))) +
      scale_colour_viridis_c("Iteration", breaks=seq(1, 10, 1)) +
      geom_line(data=observed_df, colour="orange") +
      geom_line(data=true_df, linetype=2) +
      scale_y_sqrt() +
      facet_wrap(~period) +
      ylab("Cases") +
      xlab("Onset time")
  }),
  tar_target(file_plot_cases_reporting_delay_progress, {
    
    filename <- "figures/ebola_cases_progress.pdf"
    ggsave(filename,
           plot_cases_reporting_delay_progress,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_cases_reporting_delay_progress_png, {
    
    filename <- "figures/ebola_cases_progress.png"
    ggsave(filename,
           plot_cases_reporting_delay_progress,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(plot_rt_reporting_delay_progress, {
    
    lookup <- data_sim_shortening_ebola %>% 
      select(time_onset, Rt_index) %>% 
      unique()
    
    rt_df_2 <- opt_shortening_ebola_full$Rt %>% 
      left_join(lookup,
                relationship = "many-to-many",
                by="Rt_index") %>% 
      mutate(period="E. days 1-50 & 51-90")
    rt_df_1 <- opt_shortening_ebola_full_one_value$Rt %>% 
      left_join(lookup,
                relationship = "many-to-many",
                by="Rt_index") %>% 
      mutate(period="D. days 1-90") %>% 
      filter(iteration <= 10)
    rt_df_3 <- opt_shortening_ebola_full_lots_of_vals$Rt %>% 
      left_join(lookup,
                relationship = "many-to-many",
                by="Rt_index") %>% 
      mutate(period="F. days 1-30 & 31-60 & 61-90")
    rt_df <- rt_df_1 %>% 
      bind_rows(
        rt_df_2,
        rt_df_3
      )
    
    rt_df %>% 
      filter(time_onset >= 21) %>% 
      ggplot(aes(x=time_onset, y=Rt)) +
      geom_line(aes(colour=iteration, group=as.factor(iteration))) +
      geom_line(data=resurgence_rt_df %>% filter(time_onset >= 21),
                linetype=2) +
      scale_colour_viridis_c("Iteration", breaks=seq(1, 10, 1)) +
      ylab(TeX("$R_t$")) +
      xlab("Onset time") +
      facet_wrap(~period) +
      theme(legend.position = "none",
            plot.margin = margin(r = 2.5, l=0.5, unit = "cm"))
  }),
  tar_target(plot_rt_and_cases_reporting_delay_progress, {
    plot_grid(plot_cases_reporting_delay_progress,
              plot_rt_reporting_delay_progress,
              nrow=2)
  }),
  tar_target(file_plot_rt_and_cases_reporting_delay_progress, {
    
    filename <- "figures/ebola_rt_and_cases_progress.pdf"
    ggsave(filename,
           plot_rt_and_cases_reporting_delay_progress,
           width = 10, height = 6);
    filename
  }),
  tar_target(file_plot_rt_and_cases_reporting_delay_progress_png, {
    
    filename <- "figures/ebola_rt_and_cases_progress.png"
    ggsave(filename,
           plot_rt_and_cases_reporting_delay_progress,
           width = 10, height = 6, dpi=1000);
    filename
  }),
  
  # real measles data from NL
  tar_target(measles_processed, {
    
    df <- incidenceinflation::measles_NL_2013 %>% 
      select(-date_onset) %>% 
      group_by(time_onset) %>% 
      mutate(time_reported=time_reported,
             cases_reported=cumsum(cases_reported))
    
    # filter and rebase time to start when first case appears
    df_sum <- df %>% 
      group_by(time_onset) %>% 
      summarise(cases_reported=last(cases_reported)) %>% 
      ungroup() %>% 
      filter(cases_reported>0)
    
    df %>% 
      ungroup() %>% 
      filter(time_onset >= df_sum$time_onset[1]) %>% 
      mutate(
        time_reported=time_reported - df_sum$time_onset[1] + 1,
        time_onset=time_onset - df_sum$time_onset[1] + 1
        )
  }),
  tar_target(measles_rt_index, {
    Rt_window <- 10
    days_total <- max(measles_processed$time_onset)
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
    df_Rt_index
  }),
  tar_target(measles_reporting_index, {
    
    reporting_window <- 40
    days_total <- max(measles_processed$time_onset)
    reporting_index_count <- floor(days_total / reporting_window)
    reporting_indices <- unlist(map(seq(1, reporting_index_count), ~rep(., reporting_window)))
    days_leftover <- days_total - length(reporting_indices)
    if(days_leftover > 0)
      reporting_indices <- c(reporting_indices, rep(reporting_index_count, days_leftover))
    stopifnot(length(reporting_indices)==days_total)
    
    df_reporting_index <- tibble(
      time_onset=seq(1, days_total, 1),
      reporting_piece_index=reporting_indices
    )
    df_reporting_index
  }),
  tar_target(serial_parameters_measles, {
    # from Vink et al. 2014, "Serial Intervals of Respiratory Infectious Diseases:
    # A Systematic Review and Analysis"
    # mean is overall S.I. mean found by authors
    # sd is mean of all normal sds across measles values in Table 3
    list(mean=11.7, sd=2.15)
  }),
  tar_target(results_mcmc_measles_full, 
             fit_measles(130, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, 6000)),
  tar_target(results_mcmc_measles_10, 
             fit_measles(130-10, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, niterations_measles)),
  tar_target(results_mcmc_measles_20, 
             fit_measles(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, niterations_measles)),
  tar_target(results_mcmc_measles_30, 
             fit_measles(130-30, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, 6000)),
  tar_target(results_mcmc_measles_40, 
             fit_measles(130-40, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, 6000)),
  tar_target(results_mcmc_measles_50, 
             fit_measles(130-50, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles, 6000)),
  tar_target(cases_measles_combined, {
    
    df_full <- results_mcmc_measles_full$cases %>% 
      mutate(max_time="0")
    df_10 <- results_mcmc_measles_10$cases %>% 
      mutate(max_time="10")
    df_20 <- results_mcmc_measles_20$cases %>% 
      mutate(max_time="20")
    df_30 <- results_mcmc_measles_30$cases %>% 
      mutate(max_time="30")
    df_40 <- results_mcmc_measles_40$cases %>% 
      mutate(max_time="40")
    df_50 <- results_mcmc_measles_50$cases %>% 
      mutate(max_time="50")
    
    df_full %>% 
      bind_rows(
        df_10,
        df_20,
        df_30,
        df_40,
        df_50
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "0", "10", "20", "30", "40", "50")) %>% 
      mutate(max_time=fct_rev(max_time))
  }),
  tar_target(cases_measles_true, {
    
    d <- measles_processed %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=max(cases_reported))
    
    e <- tribble(
      ~max_time,
      0,
      10,
      20,
      30,
      40,
      50) %>% 
      mutate(time_reported=130-max_time)
    
    for(i in seq_along(e$max_time)) {
      tmp <- measles_processed %>% 
        mutate(max_time=e$max_time[i]) %>% 
        filter(time_reported <= (130 - max_time))
      if(i == 1)
        big_df <- tmp
      else
        big_df <- big_df %>% bind_rows(tmp)
    }
    df_reported <- big_df %>% 
      group_by(max_time, time_onset) %>% 
      summarise(cases_reported=last(cases_reported))
    df_reported %>%
      left_join(d, by="time_onset") %>% 
      mutate(max_time=as.factor(as.character(max_time)))
  }),
  tar_target(plot_cases_measles_combined,
             plot_cases_vs_true_series_2(cases_measles_combined,
                                       cases_measles_true,
                                       measles_processed_dates,
                                       transform_type="sqrt",
                                       include_reporting=TRUE)),
  tar_target(file_plot_cases_measles_combined, {
    filename <- "figures/measles_cases.pdf"
    ggsave(filename,
           plot_cases_measles_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_cases_measles_combined_png, {
    filename <- "figures/measles_cases.png"
    ggsave(filename,
           plot_cases_measles_combined,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(summary_mcmc_measles_full,
             mcmc_summary(results_mcmc_measles_full)),
  tar_target(summary_mcmc_measles_10,
             mcmc_summary(results_mcmc_measles_10)),
  tar_target(summary_mcmc_measles_20,
             mcmc_summary(results_mcmc_measles_20)),
  tar_target(summary_mcmc_measles_30,
             mcmc_summary(results_mcmc_measles_30)),
  tar_target(summary_mcmc_measles_40,
             mcmc_summary(results_mcmc_measles_40)),
  tar_target(summary_mcmc_measles_50,
             mcmc_summary(results_mcmc_measles_50)),
  tar_target(summary_mcmc_measles_all, {
    list(
      summary_mcmc_measles_full,
      summary_mcmc_measles_10,
      summary_mcmc_measles_20,
      summary_mcmc_measles_30,
      summary_mcmc_measles_40,
      summary_mcmc_measles_50
    )
  }),
  
  # fit simple Poisson renewal model to perfect case data and compare with
  # estimates from delay method
  tar_target(model_poisson_raw, "src/stan/poisson.stan",
             format = "file"),
  tar_target(model_poisson, {
    rstan_options(auto_write=TRUE)
    stan_model(model_poisson_raw)
  }),
  tar_target(measles_generation_params, {
    
    # from Vink et al. 2014, "Serial Intervals of Respiratory Infectious Diseases:
    # A Systematic Review and Analysis"
    # mean is overall S.I. mean found by authors
    # sd is mean of all normal sds across measles values in Table 3
    mu <- 11.7
    sigma <- 2.15
    shape <- mu^2 / sigma^2
    rate <- mu / sigma^2
    
    list(
      shape=shape,
      rate=rate,
      w_max_days=40
    )
  }),
  tar_target(measles_generation_ws,
             weights_series(40, serial_parameters_measles)),
  tar_target(measles_processed_simple, {
    measles_processed %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported))
  }),
  tar_target(measles_processed_snapshot, {
    measles_processed %>% 
      filter(time_reported <= max(time_onset)) %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported))
  }),
  tar_target(fits_measles_stan_full, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple,
      measles_processed_snapshot,
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index,
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(fits_measles_stan_10, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 10)),
      measles_processed %>% 
        filter(time_reported <= (130 - 10)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-10)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(fits_measles_stan_20, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 20)),
      measles_processed %>% 
        filter(time_reported <= (130 - 20)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-20)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(fits_measles_stan_30, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 30)),
      measles_processed %>% 
        filter(time_reported <= (130 - 30)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-30)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(fits_measles_stan_40, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 40)),
      measles_processed %>% 
        filter(time_reported <= (130 - 40)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-40)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(fits_measles_stan_50, {
    fit_both_stan(
      model_poisson,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 50)),
      measles_processed %>% 
        filter(time_reported <= (130 - 50)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-50)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  
  # compare the three R estimates
  tar_target(summarise_measles_r_estimates, {
    
    # process stan fits
    R_0 <- summarise_r(fits_measles_stan_full) %>% 
      mutate(max_time=0)
    R_10 <- summarise_r(fits_measles_stan_10) %>% 
      mutate(max_time=10)
    R_20 <- summarise_r(fits_measles_stan_20) %>% 
      mutate(max_time=20)
    R_30 <- summarise_r(fits_measles_stan_30) %>% 
      mutate(max_time=30)
    R_40 <- summarise_r(fits_measles_stan_40) %>% 
      mutate(max_time=40)
    R_50 <- summarise_r(fits_measles_stan_50) %>% 
      mutate(max_time=50)
    R_stan <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_0 <- quantiles_r_other(results_mcmc_measles_full) %>% 
      mutate(max_time=0)
    R_10 <- quantiles_r_other(results_mcmc_measles_10) %>% 
      mutate(max_time=10)
    R_20 <- quantiles_r_other(results_mcmc_measles_20) %>% 
      mutate(max_time=20)
    R_30 <- quantiles_r_other(results_mcmc_measles_30) %>% 
      mutate(max_time=30)
    R_40 <- quantiles_r_other(results_mcmc_measles_40) %>% 
      mutate(max_time=40)
    R_50 <- quantiles_r_other(results_mcmc_measles_50) %>% 
      mutate(max_time=50)
    R_inflation <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(plot_measles_r_estimates, {
    
    
    df <- measles_rt_index %>% 
      left_join(summarise_measles_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      mutate(max_time=fct_relevel(max_time, "50", "40", "30", "20", "10", "0")) %>% 
      filter(time_onset >= 41)
    lookup <- measles_processed_dates %>% 
      select(time_onset, date_onset) %>% 
      unique()
    df <- df %>% 
      left_join(lookup)
    df_sum <- df %>% 
      group_by(max_time) %>% 
      summarise(
        observation_date=max(date_onset)
      )
    df <- df %>% 
      left_join(df_sum) %>% 
      mutate(observation_date=format(observation_date, "%d %b %Y")) %>% 
      mutate(observation_date=as.factor(observation_date)) %>% 
      mutate(max_time=as.numeric(as.character(max_time))) %>% 
      mutate(observation_date=fct_rev(fct_reorder(observation_date, max_time))) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    df %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
        geom_hline(yintercept=1, linetype=2) +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
        geom_line(aes(colour=variant)) +
        scale_color_manual("Estimator", values=cols) +
        scale_fill_manual("Estimator", values=cols) +
        ylab(TeX("$R_t$")) +
        xlab("Onset date") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
        facet_wrap(~observation_date)
  }),
  tar_target(file_plot_measles_r_estimates, {
    filename <- "figures/measles_rt.pdf"
    ggsave(filename,
           plot_measles_r_estimates,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_measles_r_estimates_png, {
    filename <- "figures/measles_rt.png"
    ggsave(filename,
           plot_measles_r_estimates,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(plot_measles_raw_data_delay, {
    incidenceinflation::measles_NL_2013 %>% 
      mutate(time_delay=time_reported-time_onset) %>%
      ggplot(aes(x=date_onset, y=time_delay,
                 fill=cases_reported)) +
      geom_tile() +
      scale_fill_distiller("Cases reported", direction=1,
                           breaks=c(0, 2, 4, 6, 8)) +
      xlab("Onset date") +
      ylab("Reporting delay") +
      theme(
        text=element_text(size=16)
      )
  }),
  tar_target(file_plot_measles_raw_data_delay_png, {
    filename <- "figures/measles_netherlands.png"
    ggsave(filename,
           plot_measles_raw_data_delay,
           width = 8, height = 6, dpi=1000);
    filename
  }),
  tar_target(plot_measles_delay, plot_delay_time_varying(results_mcmc_measles_full$reporting)), 
  tar_target(file_plot_measles_delay, {
    filename <- "figures/measles_delay.pdf"
    ggsave(filename,
           plot_measles_delay,
           width = 8, height = 6);
    filename
  }),
  tar_target(file_plot_measles_delay_png, {
    filename <- "figures/measles_delay.png"
    ggsave(filename,
           plot_measles_delay,
           width = 8, height = 6, dpi=1000);
    filename
  }),
  
  # simulated data with much lower case counts
  
  # simulated data with time-varying reporting period
  
  # ILI -- too few data points; I can get more but unsure what s.i. should be
  tar_target(plot_ili, {
    incidenceinflation::ILI_2014 %>% 
      mutate(time_delay=time_reported-time_onset) %>%
      ggplot(aes(x=date_onset, y=time_delay,
                 fill=cases_reported)) +
      geom_tile() +
      scale_fill_distiller("Cases reported", direction=1,
                           trans = "sqrt") +
      xlab("Onset date") +
      ylab("Reporting delay") +
      theme(
        text=element_text(size=16)
      )
  }),
  
  # dengue
  tar_target(dengue_processed, {
    
    df <- readRDS("data/raw/dengue.rds")
    a <- df %>% 
      group_by(onset_week, report_week) %>% 
      count()
    lookup <- tibble(onset_week=sort(unique(c(a$onset_week, a$report_week))),
                     time_onset=seq_along(onset_week))
    lookup1 <- lookup %>% 
      rename(time_reported=time_onset,
             report_week=onset_week)
    a <- a %>% 
      left_join(lookup, by="onset_week") %>% 
      left_join(lookup1, by="report_week")
    b <- a %>% 
      group_by(time_onset) %>% 
      mutate(cases_reported=cumsum(n))
    b
  }),
  tar_target(serial_parameters_dengue, {
    # "Estimating the effective reproduction number of dengue considering
    # temperature-dependent generation intervals" Codeco et al. 2018
    # Note this is in weeks: assume 4 weeks (the median they report)
    # with a sd of 2 (more arbitrary but would encompass their estimates)
    list(mean=4, sd=2)
  }),
  tar_target(dengue_rt_index, {
    Rt_window <- 6
    days_total <- max(dengue_processed$time_onset)
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
    df_Rt_index
  }),
  tar_target(dengue_reporting_index, {
    
    reporting_window <- 156
    days_total <- max(dengue_processed$time_onset)
    reporting_index_count <- floor(days_total / reporting_window)
    reporting_indices <- unlist(map(seq(1, reporting_index_count), ~rep(., reporting_window)))
    days_leftover <- days_total - length(reporting_indices)
    if(days_leftover > 0)
      reporting_indices <- c(reporting_indices, rep(reporting_index_count, days_leftover))
    stopifnot(length(reporting_indices)==days_total)
    
    df_reporting_index <- tibble(
      time_onset=seq(1, days_total, 1),
      reporting_piece_index=reporting_indices
    )
    df_reporting_index
  }),
  tar_target(results_mcmc_dengue_full, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= 156) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=4000)),
  tar_target(results_mcmc_dengue_10, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 10)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_20, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 20)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_30, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 30)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_40, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 40)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_50, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 50)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_60, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 60)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_70, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 70)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(results_mcmc_dengue_80, 
             fit_model(dengue_processed %>% 
                         filter(time_reported <= (156 - 80)) %>% 
                         left_join(dengue_rt_index, by="time_onset") %>% 
                         left_join(dengue_reporting_index, by="time_onset") %>% 
                         select(-c(onset_week, report_week, n)),
                       serial_parameters_dengue,
                       niterations=niterations_dengue)),
  tar_target(summary_mcmc_dengue_full,
             mcmc_summary(results_mcmc_dengue_full)),
  tar_target(summary_mcmc_dengue_10,
             mcmc_summary(results_mcmc_dengue_10)),
  tar_target(summary_mcmc_dengue_20,
             mcmc_summary(results_mcmc_dengue_20)),
  tar_target(summary_mcmc_dengue_30,
             mcmc_summary(results_mcmc_dengue_30)),
  tar_target(summary_mcmc_dengue_40,
             mcmc_summary(results_mcmc_dengue_40)),
  tar_target(summary_mcmc_dengue_50,
             mcmc_summary(results_mcmc_dengue_50)),
  tar_target(summary_mcmc_dengue_60,
             mcmc_summary(results_mcmc_dengue_60)),
  tar_target(summary_mcmc_dengue_70,
             mcmc_summary(results_mcmc_dengue_70)),
  tar_target(summary_mcmc_dengue_80,
             mcmc_summary(results_mcmc_dengue_80)),
  tar_target(summary_mcmc_dengue_all, {
    list(
      summary_mcmc_dengue_full,
      summary_mcmc_dengue_80,
      summary_mcmc_dengue_70,
      summary_mcmc_dengue_60,
      summary_mcmc_dengue_50,
      summary_mcmc_dengue_40,
      summary_mcmc_dengue_30,
      summary_mcmc_dengue_20,
      summary_mcmc_dengue_10
    )
  }),
  
  tar_target(cases_dengue_true, {
    
    d <- dengue_processed %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=max(cases_reported))
    
    e <- tribble(
      ~max_time,
      0,
      10,
      20,
      30,
      40,
      50,
      60,
      70,
      80) %>% 
      mutate(time_reported=156-max_time)
    
    for(i in seq_along(e$max_time)) {
      tmp <- dengue_processed %>% 
        mutate(max_time=e$max_time[i]) %>% 
        filter(time_reported <= (156 - max_time))
      if(i == 1)
        big_df <- tmp
      else
        big_df <- big_df %>% bind_rows(tmp)
    }
    df_reported <- big_df %>% 
      group_by(max_time, time_onset) %>% 
      summarise(cases_reported=last(cases_reported))
    df_reported %>%
      left_join(d, by="time_onset") %>% 
      mutate(max_time=as.factor(as.character(max_time)))
  }),
  tar_target(cases_dengue_combined, {
    
    df_full <- results_mcmc_dengue_full$cases %>% 
      mutate(max_time="0")
    df_10 <- results_mcmc_dengue_10$cases %>% 
      mutate(max_time="10")
    df_20 <- results_mcmc_dengue_20$cases %>% 
      mutate(max_time="20")
    df_30 <- results_mcmc_dengue_30$cases %>% 
      mutate(max_time="30")
    df_40 <- results_mcmc_dengue_40$cases %>% 
      mutate(max_time="40")
    df_50 <- results_mcmc_dengue_50$cases %>% 
      mutate(max_time="50")
    df_60 <- results_mcmc_dengue_60$cases %>% 
      mutate(max_time="60")
    df_70 <- results_mcmc_dengue_70$cases %>% 
      mutate(max_time="70")
    df_80 <- results_mcmc_dengue_80$cases %>% 
      mutate(max_time="80")
    
    df_full %>% 
      bind_rows(
        df_10,
        df_20,
        df_30,
        df_40,
        df_50,
        df_60,
        df_70,
        df_80
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "0", "10", "20", "30", "40", "50", "60", "70", "80")) %>% 
      mutate(max_time=fct_rev(max_time))
  }),
  tar_target(plot_cases_dengue_combined,
             plot_cases_vs_true_series_one_line(cases_dengue_combined,
                                       cases_dengue_true,
                                       dengue_processed
                                       )),
  tar_target(file_plot_cases_dengue_combined, {
    filename <- "figures/dengue_cases.pdf"
    ggsave(filename,
           plot_cases_dengue_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_cases_dengue_combined_png, {
    filename <- "figures/dengue_cases.png"
    ggsave(filename,
           plot_cases_dengue_combined,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(dengue_processed_simple, {
    dengue_processed %>% 
      group_by(time_onset) %>% 
      summarise(C=last(cases_reported))
  }),
  tar_target(dengue_generation_params, {
    # "Estimating the effective reproduction number of dengue considering
    # temperature-dependent generation intervals" Codeco et al. 2018
    # Note this is in weeks: assume 4 weeks (the median they report)
    # with a sd of 2 (more arbitrary but would encompass their estimates)
    mu <- 4
    sigma <- 2
    shape <- mu^2 / sigma^2
    rate <- mu / sigma^2
    
    list(
      shape=shape,
      rate=rate,
      w_max_days=40
    )
  }),
  tar_target(dengue_generation_ws,
             weights_series(40, serial_parameters_dengue)),
  tar_target(fits_dengue_stan_0, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156, 2000, model_poisson)),
  tar_target(fits_dengue_stan_10, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-10, 2000, model_poisson)),
  tar_target(fits_dengue_stan_20, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-20, 2000, model_poisson)),
  tar_target(fits_dengue_stan_30, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-30, 2000, model_poisson)),
  tar_target(fits_dengue_stan_40, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-40, 2000, model_poisson)),
  tar_target(fits_dengue_stan_50, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-50, 2000, model_poisson)),
  tar_target(fits_dengue_stan_60, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-60, 2000, model_poisson)),
  tar_target(fits_dengue_stan_70, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-70, 2000, model_poisson)),
  tar_target(fits_dengue_stan_80, fit_both_stan_dengue(
    dengue_processed_simple, dengue_processed, dengue_generation_ws, dengue_rt_index, 156-80, 2000, model_poisson)),
  tar_target(summarise_dengue_r_estimates, {
    
    # process stan fits
    R_0 <- summarise_r(fits_dengue_stan_0) %>% 
      mutate(max_time=0)
    R_10 <- summarise_r(fits_dengue_stan_10) %>% 
      mutate(max_time=10)
    R_20 <- summarise_r(fits_dengue_stan_20) %>% 
      mutate(max_time=20)
    R_30 <- summarise_r(fits_dengue_stan_30) %>% 
      mutate(max_time=30)
    R_40 <- summarise_r(fits_dengue_stan_40) %>% 
      mutate(max_time=40)
    R_50 <- summarise_r(fits_dengue_stan_50) %>% 
      mutate(max_time=50)
    R_60 <- summarise_r(fits_dengue_stan_60) %>% 
      mutate(max_time=60)
    R_70 <- summarise_r(fits_dengue_stan_70) %>% 
      mutate(max_time=70)
    R_80 <- summarise_r(fits_dengue_stan_80) %>% 
      mutate(max_time=80)
    R_stan <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50,
        R_60,
        R_70,
        R_80
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_0 <- quantiles_r_other(results_mcmc_dengue_full) %>% 
      mutate(max_time=0)
    R_10 <- quantiles_r_other(results_mcmc_dengue_10) %>% 
      mutate(max_time=10)
    R_20 <- quantiles_r_other(results_mcmc_dengue_20) %>% 
      mutate(max_time=20)
    R_30 <- quantiles_r_other(results_mcmc_dengue_30) %>% 
      mutate(max_time=30)
    R_40 <- quantiles_r_other(results_mcmc_dengue_40) %>% 
      mutate(max_time=40)
    R_50 <- quantiles_r_other(results_mcmc_dengue_50) %>% 
      mutate(max_time=50)
    R_60 <- quantiles_r_other(results_mcmc_dengue_60) %>% 
      mutate(max_time=60)
    R_70 <- quantiles_r_other(results_mcmc_dengue_70) %>% 
      mutate(max_time=70)
    R_80 <- quantiles_r_other(results_mcmc_dengue_80) %>% 
      mutate(max_time=80)
    R_inflation <- R_0 %>% 
      bind_rows(
        R_10,
        R_20,
        R_30,
        R_40,
        R_50,
        R_60,
        R_70,
        R_80
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(plot_dengue_r_estimates, {
    
    lookup <- dengue_processed %>% 
      select(time_onset, onset_week) %>% 
      unique()
    
    # necessary since Rt indices extend beyond obs date occassionally whereas cases don't
    test_df <- cases_dengue_combined %>% 
      select(time_onset, max_time) %>% 
      unique() %>% 
      mutate(is_after=1)
    
    # complications from same issue above
    lookup_observation_time <- dengue_rt_index %>% 
      left_join(summarise_dengue_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      filter(time_onset <= 156) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      left_join(lookup) %>% 
      right_join(test_df) %>% 
      arrange(desc(max_time), time_onset) %>% 
      group_by(max_time) %>% 
      mutate(observation_point=last(onset_week)) %>% 
      mutate(observation_point=format(observation_point, "%d %b %Y")) %>% 
      mutate(observation_point=as.factor(observation_point)) %>% 
      mutate(observation_point=fct_rev(observation_point)) %>% 
      ungroup() %>% 
      select(max_time, observation_point) %>% 
      unique()
    
    df <- dengue_rt_index %>% 
      left_join(summarise_dengue_r_estimates,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      filter(time_onset <= 156) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      left_join(lookup) %>% 
      left_join(lookup_observation_time) %>% 
      filter(time_onset >= 60) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    ggplot(df, aes(x=onset_week, y=middle, fill=variant)) +
      geom_hline(yintercept=1, linetype=2) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~observation_point) +
      scale_y_sqrt() +
      scale_x_date(date_breaks = "6 month", date_labels = "%b %Y")
  }),
  tar_target(file_plot_dengue_r_estimates, {
    filename <- "figures/dengue_rt.pdf"
    ggsave(filename,
           plot_dengue_r_estimates,
           width = 10, height = 6);
    filename
  }),
  tar_target(file_plot_dengue_r_estimates_png, {
    filename <- "figures/dengue_rt.png"
    ggsave(filename,
           plot_dengue_r_estimates,
           width = 10, height = 6, dpi=1000);
    filename
  }),
  
  
  
  # using stan version
  tar_target(stan_model, { # not currently allowing me to pass this in: some sort of Rcpp thing
    cmdstanr::cmdstan_model("src/stan/conditional_renewal.stan",
                            compile_model_methods=TRUE,
                            force_recompile=TRUE)
  }),
  tar_target(opt_mcmc_measles_20, 
             fit_measles_uncertain_cases_stan(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                         serial_parameters_measles,
                         15,
                         is_negative_binomial=FALSE,
                         is_maximise=TRUE)),
  tar_target(fit_mcmc_measles_20,
             {
               initial_cases_true <- opt_mcmc_measles_20$cases %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_Rt <- opt_mcmc_measles_20$Rt %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_reporting_parameters <- opt_mcmc_measles_20$reporting %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)

               fit_measles_uncertain_cases_stan(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                                              serial_parameters_measles,
                                              600,
                                              is_negative_binomial=FALSE,
                                              is_maximise=FALSE,
                                              initial_cases_true=initial_cases_true,
                                              initial_reporting_parameters=initial_reporting_parameters,
                                              initial_Rt=initial_Rt)
             }),
  tar_target(fit_mcmc_measles_20_short,
             {
               initial_cases_true <- opt_mcmc_measles_20$cases %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_Rt <- opt_mcmc_measles_20$Rt %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               print(initial_Rt)
               initial_reporting_parameters <- opt_mcmc_measles_20$reporting %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               
               fit_measles_uncertain_cases_stan(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                                                serial_parameters_measles,
                                                4000,
                                                is_negative_binomial=FALSE,
                                                is_maximise=FALSE,
                                                initial_cases_true=initial_cases_true,
                                                initial_reporting_parameters=initial_reporting_parameters,
                                                initial_Rt=initial_Rt)
             }),
  tar_target(opt_mcmc_measles_20_nb, 
             fit_measles_uncertain_cases_stan(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                                              serial_parameters_measles,
                                              15,
                                              is_negative_binomial=TRUE,
                                              is_maximise=TRUE)),
  tar_target(fit_mcmc_measles_20_short_nb,
             {
               initial_cases_true <- opt_mcmc_measles_20_nb$cases %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_Rt <- opt_mcmc_measles_20_nb$Rt %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_reporting_parameters <- opt_mcmc_measles_20_nb$reporting %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_overdispersion <- opt_mcmc_measles_20_nb$other %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration) %>% 
                 pull(overdispersion)
               
               fit_measles_uncertain_cases_stan(130-20, measles_processed, measles_rt_index, measles_reporting_index,
                                                serial_parameters_measles,
                                                100,
                                                is_negative_binomial=TRUE,
                                                is_maximise=FALSE,
                                                initial_cases_true=initial_cases_true,
                                                initial_reporting_parameters=initial_reporting_parameters,
                                                initial_Rt=initial_Rt,
                                                initial_overdispersion=initial_overdispersion,
                                                use_stan_sampling=TRUE,
                                                n_stan_warmup=50,
                                                n_stan_iterations=10,
                                                step_size=1e-3)
             }),
  tar_target(summary_mcmc_measles_20_short_nb,
             mcmc_summary_stan_uncertain_cases(fit_mcmc_measles_20_short_nb)),
  tar_target(opt_uncertain_cases_stan_step_down_up,
             fit_model_uncertain_cases_stan(data_sim_step_down_up %>% 
                                    select(-cases_true) %>% 
                                    filter(time_reported <=50),
                                  serial_parameters_covid,
                                  niterations=10,
                                  is_maximise=TRUE)),
  tar_target(fit_uncertain_cases_stan_step_down_up, # this looks ok as approximation
             {
               
             initial_cases_true <- opt_uncertain_cases_stan_step_down_up$cases %>% 
               filter(iteration==max(iteration)) %>% 
               select(-iteration)
             initial_Rt <- opt_uncertain_cases_stan_step_down_up$Rt %>% 
               filter(iteration==max(iteration)) %>% 
               select(-iteration)
             initial_reporting_parameters <- opt_uncertain_cases_stan_step_down_up$reporting %>% 
               filter(iteration==max(iteration)) %>% 
               select(-iteration)
             fit_model_uncertain_cases_stan(data_sim_step_down_up %>% 
                                              select(-cases_true) %>% 
                                              filter(time_reported <=50),
                                            serial_parameters_covid,
                                            niterations=200,
                                            is_maximise=FALSE,
                                            initial_cases_true=initial_cases_true,
                                            initial_reporting_parameters=initial_reporting_parameters,
                                            initial_Rt=initial_Rt)
               }),
  tar_target(opt_resurgence_full_nb,
             fit_ebola_uncertain_cases_stan(125,
                                            data_sim_resurgence_ebola,
                                            serial_parameters_ebola,
                                            10,
                                            stan_model,
                                            is_negative_binomial=TRUE,
                                            is_maximise=TRUE)),
  tar_target(fit_resurgence_full_nb,
             {
               
               initial_cases_true <- opt_resurgence_full_nb$cases %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_Rt <- opt_resurgence_full_nb$Rt %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_reporting_parameters <- opt_resurgence_full_nb$reporting %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration)
               initial_overdispersion <- opt_resurgence_full_nb$other %>% 
                 filter(iteration==max(iteration)) %>% 
                 select(-iteration) %>% 
                 pull(overdispersion)
               fit_ebola_uncertain_cases_stan(125,
                                              data_sim_resurgence_ebola,
                                              serial_parameters_ebola,
                                              600,
                                              stan_model,
                                              is_negative_binomial=TRUE,
                                              is_maximise=FALSE,
                                              initial_cases_true=initial_cases_true,
                                              initial_reporting_parameters=initial_reporting_parameters,
                                              initial_Rt=initial_Rt,
                                              initial_overdispersion=initial_overdispersion)
             }),
  tar_target(fit_resurgence_full_po, 
             fit_model(data_sim_resurgence_ebola %>% 
                         select(-cases_true) %>%
                         filter(time_reported <= 125),
                       serial_parameters_ebola,
                       niterations=200,
                       is_negative_binomial=FALSE)),
  tar_target(fits_measles_stan_50_nb, {
    fit_both_stan(
      model_nb,
      measles_processed_simple %>% 
        filter(time_onset <= (130 - 50)),
      measles_processed %>% 
        filter(time_reported <= (130 - 50)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      measles_generation_ws,
      measles_rt_index$Rt_index[1:(130-50)],
      I_init=rep(1, 40),
      n_iterations=200)
  }),
  tar_target(opt_measles_50_nb, {
    fit_model_uncertain_cases_stan(
      measles_processed %>% 
        filter(time_reported <= (130 - 50)) %>% 
        thin_series() %>% 
        left_join(measles_rt_index, by="time_onset") %>% 
        left_join(measles_reporting_index, by="time_onset"),
      serial_parameters_measles,
      niterations=20,
      is_negative_binomial=TRUE,
      is_maximise=TRUE)
  }),
  tar_target(fit_measles_50_nb, {
    
    initial_cases_true <- opt_measles_50_nb$cases %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_Rt <- opt_measles_50_nb$Rt %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_reporting_parameters <- opt_measles_50_nb$reporting %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration)
    initial_overdispersion <- opt_measles_50_nb$other %>% 
      filter(iteration==max(iteration)) %>% 
      select(-iteration) %>% 
      pull(overdispersion)
    
    fit_model_uncertain_cases_stan(
      measles_processed %>% 
        filter(time_reported <= (130 - 50)) %>% 
        thin_series() %>% 
        left_join(measles_rt_index, by="time_onset") %>% 
        left_join(measles_reporting_index, by="time_onset"),
      serial_parameters_measles,
      niterations=200,
      is_negative_binomial=TRUE,
      is_maximise=FALSE,
      initial_cases_true=initial_cases_true,
      initial_reporting_parameters=initial_reporting_parameters,
      initial_Rt=initial_Rt,
      initial_overdispersion=initial_overdispersion,
      use_stan_sampling=TRUE,
      step_size=1e-3)
    
  }),
  tar_target(opt_measles_40_nb, {
    fit_model_uncertain_cases_stan(
      measles_processed %>% 
        filter(time_reported <= (130 - 40)) %>% 
        thin_series() %>% 
        left_join(measles_rt_index, by="time_onset") %>% 
        left_join(measles_reporting_index, by="time_onset"),
      serial_parameters_measles,
      niterations=20,
      is_negative_binomial=TRUE,
      is_maximise=TRUE)
  }),
  tar_target(plot_measles_nb_vs_po, {
    
    cases_nb <- fit_mcmc_measles_20_short_nb$cases
    max_iteration <- max(cases_nb$iteration)
    cases_nb <- cases_nb %>% 
      filter(iteration >= max_iteration / 2) %>% 
      mutate(model="negative binomial")
    cases_po <- results_mcmc_measles_20$cases
    max_iteration <- max(cases_po$iteration)
    cases_po <- cases_po %>% 
      filter(iteration >= max_iteration / 2) %>% 
      mutate(model="Poisson")
    cases_both <- cases_po %>% 
      bind_rows(cases_nb)
    cases_sum <- cases_both %>% 
      group_by(model, time_onset) %>% 
      summarise(
        lower=quantile(cases_true, 0.025),
        middle=quantile(cases_true, 0.5),
        upper=quantile(cases_true, 0.975)
      )
    
    lookup <- measles_processed_dates %>% 
      select(date_onset, time_onset) %>% 
      unique()
    cases_sum <- cases_sum %>% 
      left_join(lookup)
    
    max_date <- max(cases_sum$date_onset)
    df_true <- measles_processed_dates %>% 
      group_by(date_onset) %>% 
      summarise(middle=last(cases_reported)) %>% 
      ungroup() %>% 
      filter(date_onset<=max_date)
    
    cases_sum %>% 
      ggplot(aes(x=date_onset, y=middle)) +
      geom_line(aes(colour=model, group=as.factor(model))) +
      geom_ribbon(aes(ymin=lower, ymax=upper, fill=model, group=as.factor(model)), alpha=0.4) +
      geom_line(data=df_true) +
      scale_color_brewer("Model", palette = "Dark2") +
      scale_fill_brewer("Model", palette = "Dark2") +
      xlab("Onset date") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      ylab("Cases")
  }),
  tar_target(file_plot_measles_nb_vs_po, {
    filename <- "figures/meales_nb_vs_poisson.pdf"
    ggsave(filename,
           plot_measles_nb_vs_po,
           width = 8, height = 4);
    filename
  }),
  tar_target(file_plot_measles_nb_vs_po_png, {
    filename <- "figures/meales_nb_vs_poisson.png"
    ggsave(filename,
           plot_measles_nb_vs_po,
           width = 8, height = 4, dpi=1000);
    filename
  }),
  
  # A waning epidemic
  tar_target(rt_waning, {
    v_Rt <- c(rep(1.2, 40), rep(0.5, 60))
    v_Rt
  }),
  tar_target(rt_waning_df, {
    tibble(Rt=rt_waning,
           time_onset=seq_along(Rt))
  }),
  tar_target(data_sim_waning,
             generate_data(rt_waning,
                           serial_parameters_ebola,
                           reporting_parameters_10_3)),
  tar_target(prob_waning_over_gold_standard, {
    
    cases_recent <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=first(cases_true)) %>% 
      pull(cases_true)
    cases_recent <- rev(cases_recent)
    cases_recent <- cases_recent[1:length(ebola_generation_ws)]
    
    prob_outbreak_over(50, 0.5, cases_recent, ebola_generation_ws)
    }),
  tar_target(prob_waning_over_naive, {
    
    cases_recent <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=last(cases_reported)) %>% 
      pull(cases_true)
    cases_recent <- rev(cases_recent)
    cases_recent <- cases_recent[1:length(ebola_generation_ws)]
    
    prob_outbreak_over(50, 0.5, cases_recent, ebola_generation_ws)
  }),
  tar_target(fit_waning_mcmc,
             fit_model(data_sim_waning %>%
                         mutate(reporting_piece_index=1) %>% 
                         select(-cases_true),
                       serial_parameters_ebola,
                       niterations=200)),
  tar_target(prob_waning_over_incidenceinflation,
             prob_outbreak_over_incidence_inflation(
               50, fit_waning_mcmc, ebola_generation_ws,
               n_iterates=1000)),
  tar_target(projections_waning_over_gold_standard, {
    
    cases_recent <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=first(cases_true)) %>% 
      pull(cases_true)
    cases_recent <- rev(cases_recent)
    cases_recent <- cases_recent[1:length(ebola_generation_ws)]
    
    forward_simulations_many(50, 0.5, cases_recent, ebola_generation_ws,
                             max(data_sim_waning$time_onset))
  }),
  tar_target(projections_waning_over_naive, {
    
    cases_recent <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(cases_true=last(cases_reported)) %>% 
      pull(cases_true)
    cases_recent <- rev(cases_recent)
    cases_recent <- cases_recent[1:length(ebola_generation_ws)]
    
    forward_simulations_many(50, 0.5, cases_recent, ebola_generation_ws,
                             max(data_sim_waning$time_onset))
  }),
  tar_target(projections_waning_over_incidenceinflation, {
    forward_simulations_many_incidenceinflation(50, fit_waning_mcmc, ebola_generation_ws,
                                                max(data_sim_waning$time_onset))
  }),
  tar_target(plot_projections_waning_over, {
    
    df_gold_standard <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(I=first(cases_true)) %>% 
      expand_grid(iteration=1:max(projections_waning_over_gold_standard$iteration))
    df_naive <- data_sim_waning %>% 
      group_by(time_onset) %>% 
      summarise(I=last(cases_reported)) %>% 
      expand_grid(iteration=1:max(projections_waning_over_gold_standard$iteration))
    cases_df <- fit_waning_mcmc$cases
    max_iteration <- max(cases_df$iteration)
    df_fit_mcmc <- cases_df %>% 
      filter(iteration >= max_iteration / 2) %>% 
      mutate(combined=paste0(iteration, chain)) %>% 
      mutate(counter=as.numeric(as.factor(combined)))
    max_projections_iteration <- max(projections_waning_over_incidenceinflation$iteration)
    for(i in 1:max_projections_iteration) {
      id <- sample(1:max(df_fit_mcmc$counter), 1)
      df_tmp <- df_fit_mcmc %>%
        filter(counter==id) %>% 
        rename(I=cases_true) %>% 
        select(time_onset, I) %>% 
        mutate(iteration=i)
      if(i == 1)
        df_incidenceinflation <- df_tmp
      else
        df_incidenceinflation <- df_incidenceinflation %>% bind_rows(df_tmp)
    }
     
    df_1 <- df_gold_standard %>% 
      mutate(case_type=if_else(time_onset==100, "future", "past")) %>% 
      bind_rows(projections_waning_over_gold_standard %>% mutate(case_type="future")) %>% 
      mutate(type="A. gold standard") %>% 
      arrange(iteration, time_onset)
    df_2 <- df_naive %>% 
      mutate(case_type=if_else(time_onset==100, "future", "past")) %>% 
      bind_rows(projections_waning_over_naive %>% mutate(case_type="future")) %>% 
      mutate(type="B. naive") %>% 
      arrange(iteration, time_onset)
    df_3 <- df_incidenceinflation %>%
      mutate(case_type=if_else(time_onset==100, "future", "past")) %>% 
      bind_rows(projections_waning_over_incidenceinflation %>% mutate(case_type="future")) %>% 
      mutate(type="C. incidenceinflation") %>% 
      arrange(iteration, time_onset)
    
    cols <- c("past" = "black", "future"="blue")
    df_1 %>% 
      bind_rows(
        df_2,
        df_3
      ) %>% 
      filter(time_onset >= 70) %>% 
      mutate(type=as.factor(type)) %>% 
      mutate(type=fct_relevel(type, "A. gold standard", "B. naive", "C. incidenceinflation")) %>% 
      ggplot(aes(x=time_onset, y=I, group=as.factor(iteration))) +
      geom_line(alpha=0.3, aes(colour=case_type)) +
      facet_wrap(~type) +
      xlab("Onset time") +
      ylab("Cases") +
      scale_colour_manual(values=cols) +
      geom_vline(xintercept = 100, linetype=2) +
      theme(legend.position = "none")
  }),
  tar_target(file_plot_projections_waning_over, {
    filename <- "figures/projections_waning.pdf"
    ggsave(filename,
           plot_projections_waning_over,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_projections_waning_over_png, {
    filename <- "figures/projections_waning.png"
    ggsave(filename,
           plot_projections_waning_over,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  
  ## examining measles data in more depth
  tar_target(plot_measles_ecdf, {
    
    df_summary <- measles_processed %>% 
      group_by(time_onset) %>% 
      filter(max(cases_reported) > 5) %>% 
      mutate(
        ecdf=cases_reported/max(cases_reported),
        reporting_delay=seq_along(time_onset)
      )
    df_sum <- df_summary %>% 
      mutate(onset_period=case_when(
        time_onset<50~"(i) <50 days",
        time_onset>=50 & time_onset<100~"(ii) >=50 days & < 100 days",
        time_onset>=100~"(iii) >= 100 days"
      )) %>% 
      mutate(onset_period=as.factor(onset_period)) %>% 
      mutate(onset_period=fct_relevel(onset_period, "(i) <50 days", "(ii) >=50 days & < 100 days", "(iii) >= 100 days")) %>% 
      group_by(time_onset) %>% 
      mutate(max_cases_reported=max(cases_reported)) %>% 
      mutate(cases_observed=c(0, diff(cases_reported))) %>% 
      ungroup()
    df_overall <- df_sum %>% 
      group_by(onset_period, reporting_delay) %>% 
      summarise(ecdf=mean(ecdf))
    onset_periods <- unique(df_overall$onset_period)
    d50s <- vector(length = length(onset_periods))
    for(i in seq_along(onset_periods)) {
      df_tmp <- df_overall %>% filter(onset_period==onset_periods[i])
      f_iecdf <- approxfun(df_tmp$ecdf, df_tmp$reporting_delay)
      d50s[i] <- f_iecdf(0.5)
    }
    
    df_50s <- tibble(reporting_delay=d50s,
                     ecdf=0.5,
                     onset_period=onset_periods)
      
    df_sum %>% 
      ggplot(aes(x=reporting_delay, y=ecdf)) +
      geom_line(alpha=0.4, aes(group=as.factor(time_onset),
                               colour=max_cases_reported)) +
      geom_line(data=df_overall, colour="orange", linetype=2) +
      geom_segment(data=df_50s, aes(xend=1, yend=ecdf),
                   colour="orange") +
      geom_segment(data=df_50s, aes(xend=reporting_delay, yend=0),
                   colour="orange") +
      facet_wrap(~onset_period) +
      scale_x_sqrt() +
      scale_colour_viridis_c("Max\ncases") +
      xlab("Reporting delay, days") +
      ylab("eCDF")
  }),
  tar_target(file_plot_measles_ecdf, {
    filename <- "figures/measles_ecdf.pdf"
    ggsave(filename,
           plot_measles_ecdf,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_measles_ecdf_png, {
    filename <- "figures/measles_ecdf.png"
    ggsave(filename,
           plot_measles_ecdf,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(measles_processed_dates, {
    df <- incidenceinflation::measles_NL_2013 %>% 
      group_by(time_onset) %>% 
      mutate(time_reported=time_reported,
             cases_reported=cumsum(cases_reported))
    
    # filter and rebase time to start when first case appears
    df_sum <- df %>% 
      group_by(time_onset) %>% 
      summarise(cases_reported=last(cases_reported)) %>% 
      ungroup() %>% 
      filter(cases_reported>0)
    
    df <- df %>% 
      ungroup() %>% 
      filter(time_onset >= df_sum$time_onset[1]) %>% 
      mutate(
        time_reported=time_reported - df_sum$time_onset[1] + 1,
        time_onset=time_onset - df_sum$time_onset[1] + 1
      )
  }),
  tar_target(plot_measles_cases_simple, {
    
    df <- measles_processed_dates
    
    lookup <- df %>% 
      select(time_onset, date_onset) %>% 
      unique()
    
    df_rect <- tribble(
      ~xmin, ~xmax, ~period,
      1, 50, 1,
      50, 100, 2,
      100, max(df$time_onset), 3
    ) %>% 
      pivot_longer(-period) %>% 
      rename(time_onset=value) %>% 
      left_join(lookup) %>% 
      select(-time_onset) %>% 
      pivot_wider(names_from = name,
                  values_from = date_onset) %>% 
      mutate(ymin=0,
             ymax=max(df$cases_reported))
    df_text <- tribble(
      ~time_onset, ~label,
      25, "(i)",
      75, "(ii)",
      115, "(iii)"
    ) %>%
      left_join(lookup) %>% 
      select(-time_onset) %>% 
      mutate(y=35)
    
    df %>% 
      group_by(time_onset, date_onset) %>% 
      summarise(cases_true=last(cases_reported)) %>% 
      ggplot() +
      geom_rect(data=df_rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                  fill=as.factor(period)),
                alpha=0.25) +
      geom_line(aes(x=date_onset, y=cases_true)) +
      geom_text(data=df_text, aes(x=date_onset, y=y, label=label)) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      xlab("Onset date") +
      ylab("Cases") +
      scale_fill_brewer(palette = "Dark2") +
      theme(
        legend.position = "none"
      )
  }),
  tar_target(plot_measles_cases_ecdf, {
    plot_grid(plot_measles_cases_simple,
              plot_measles_ecdf,
              labels=c("A.", "B."),
              nrow=2)
  }),
  tar_target(file_plot_measles_cases_ecdf, {
    filename <- "figures/measles_cases_ecdf.pdf"
    ggsave(filename,
           plot_measles_cases_ecdf,
           width = 10, height = 8);
    filename
  }),
  tar_target(file_plot_measles_cases_ecdf_png, {
    filename <- "figures/measles_cases_ecdf.png"
    ggsave(filename,
           plot_measles_cases_ecdf,
           width = 10, height = 8, dpi=1000);
    filename
  }),
  
  tar_target(stan_delay_raw, "src/stan/delay_period_distribution.stan",
             format = "file"),
  tar_target(stan_model_delay, {
    rstan_options(auto_write=TRUE)
    rstan::stan_model(stan_delay_raw)
    }),
  tar_target(measles_delay_stan, {
    
    df_tmp <- measles_processed %>% 
      incidenceinflation::thin_series() %>% 
      group_by(time_onset) %>% 
      mutate(cases_true=max(cases_reported))
    
    reporting_days_sim <- seq(0, 40, 0.1)
    
    data_stan <- list(
      N_delay=nrow(df_tmp),
      time_reported=df_tmp$time_reported,
      time_onset=df_tmp$time_onset,
      cases_reported=df_tmp$cases_reported,
      cases_true=df_tmp$cases_true,
      n_reporting_window=1,
      reporting_window=rep(1, length(df_tmp$time_onset)),
      is_gamma_reporting_delay=0,
      reporting_days_sim=reporting_days_sim,
      N_sim=length(reporting_days_sim),
      is_binomial=1
    )
    
    options(mc.cores=4)
    init_fn <- function() {
      list(
        theta=matrix(c(10, 5), nrow = 1)
        )
    }
    optimizing(stan_model_delay, data=data_stan, init=init_fn, as_vector=FALSE)
  }),
  tar_target(measles_delay_stan_beta_binomial, {
    
    df_tmp <- measles_processed %>% 
      incidenceinflation::thin_series() %>% 
      group_by(time_onset) %>% 
      mutate(cases_true=max(cases_reported))
    
    reporting_days_sim <- seq(0, 40, 0.1)
    
    data_stan <- list(
      N_delay=nrow(df_tmp),
      time_reported=df_tmp$time_reported,
      time_onset=df_tmp$time_onset,
      cases_reported=df_tmp$cases_reported,
      cases_true=df_tmp$cases_true,
      n_reporting_window=1,
      reporting_window=rep(1, length(df_tmp$time_onset)),
      is_gamma_reporting_delay=0,
      reporting_days_sim=reporting_days_sim,
      N_sim=length(reporting_days_sim),
      is_binomial=0
    )
    
    options(mc.cores=4)
    init_fn <- function() {
      list(
        theta=matrix(c(10, 5), nrow = 1),
        kappa=as.array(c(6))
      )
    }
    optimizing(stan_model_delay, data=data_stan, init=init_fn, as_vector=FALSE)
  }),
  tar_target(plot_measles_delay_stan_ecdf, {
    
    df_tmp <- measles_processed %>% 
      incidenceinflation::thin_series() %>% 
      group_by(time_onset) %>% 
      mutate(cases_true=max(cases_reported))
    
    fit <- measles_delay_stan
    v_ecdf_binomial <- fit$par$v_ecdf
    fit <- measles_delay_stan_beta_binomial
    v_ecdf_beta_binomial <- fit$par$v_ecdf
    v_cases_remaining <- fit$par$v_cases_remaining
    v_cases_observed <- fit$par$v_cases_observed
    
    df_short <- tibble(
      time_onset=df_tmp$time_onset[1:(nrow(df_tmp) - 1)],
      time_reported=df_tmp$time_reported[1:(nrow(df_tmp) - 1)],
      cases_true=df_tmp$cases_true[1:(nrow(df_tmp) - 1)],
      cases_remaining=v_cases_remaining,
      cases_observed=v_cases_observed,
      binomial=v_ecdf_binomial,
      beta_binomial=v_ecdf_beta_binomial
    ) %>% 
      group_by(time_onset) %>% 
      mutate(reporting_delay=seq_along(time_onset)) %>% 
      filter(cases_remaining > 0) %>% 
      filter(cases_remaining != cases_observed) %>% # justified because each onset time always ends with one of these and we've filtered out the zero observations
      mutate(fraction_observed=cases_observed/cases_remaining)
    
    df_longer <- df_short %>% 
      pivot_longer(c(beta_binomial, binomial)) %>% 
      mutate(name=if_else(name=="beta_binomial", "B. beta-binomial", "A. binomial")) %>% 
      mutate(name=as.factor(name)) %>% 
      mutate(name=fct_relevel(name, "A. binomial", "B. beta-binomial"))
    
    ggplot(df_longer, aes(x=value)) +
      geom_histogram() +
      xlab("eCDF") +
      ylab("Count") +
      facet_wrap(~name) +
      theme(
        text=element_text(size=16)
      )
  }),
  tar_target(file_plot_measles_ecdfs_binomial_beta, {
    filename <- "figures/measles_beta_binomial.pdf"
    ggsave(filename,
           plot_measles_delay_stan_ecdf,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_measles_ecdfs_binomial_beta_png, {
    filename <- "figures/measles_beta_binomial.png"
    ggsave(filename,
           plot_measles_delay_stan_ecdf,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(measles_largest_cases, {
    
    df_big_cases <- measles_processed %>%
      filter(time_onset >= 20) %>% 
      group_by(time_onset) %>% 
      summarise(cases_reported=max(cases_reported)) %>% 
      mutate(is_big=if_else(cases_reported > 35, TRUE, FALSE)) %>% # corresponds to day with most cases
      ungroup() %>% 
      select(-cases_reported)
    
    df_big <- measles_processed %>%
      filter(time_onset >= 20) %>% 
      left_join(df_big_cases) %>% 
      filter(is_big) %>% 
      group_by(time_onset) %>% 
      mutate(time_reported=seq_along(time_reported))
    
    df_big
  }),
  tar_target(fit_measles_largest, {
    
    reporting_days_sim <- seq(0, 40, 0.1)
    
    df_tmp <- measles_largest_cases
    data_stan <- list(
      N_delay=nrow(df_tmp),
      time_reported=df_tmp$time_reported,
      time_onset=rep(1, nrow(df_tmp)),
      cases_reported=df_tmp$cases_reported,
      cases_true=rep(max(df_tmp$cases_reported), nrow(df_tmp)),
      n_reporting_window=1,
      reporting_window=rep(1, length(df_tmp$time_onset)),
      is_gamma_reporting_delay=1,
      reporting_days_sim=reporting_days_sim,
      N_sim=length(reporting_days_sim),
      is_binomial=1
    )
    
    options(mc.cores=4)
    init_fn <- function() {
      list(
        theta=matrix(c(10, 5), nrow = 1)
      )
    }
    options(mc.cores=4)
    sampling(stan_model_delay, data=data_stan, init=init_fn, iter=200, chains=4)
  }),
  tar_target(measles_largest_sim_cases, {
    
    fit <- fit_measles_largest
    cases_reported <- rstan::extract(fit, "v_cases_reported")[[1]] %>% 
      as.data.frame() %>% 
      mutate(iteration=seq_along(V1)) %>% 
      pivot_longer(-iteration) %>% 
      group_by(iteration) %>% 
      mutate(time_reported=seq_along(iteration)) %>% 
      rename(cases_reported=value) %>% 
      ungroup() %>% 
      select(-name)
    
    cases_reported
  }),
  tar_target(plot_measles_ecdf_absolute, {
    
    measles_largest_cases %>%   
      ggplot(aes(x=time_reported,
                 y=cases_reported)) +
      geom_line(data=measles_largest_sim_cases,
                aes(group=as.factor(iteration))) +
      geom_line(colour="orange") +
      xlab("Reporting delay, days") +
      ylab("Cases reported") +
      theme(legend.position = c(0.7, 0.75))
  }),
  tar_target(measles_largest_sim_cases_observed, {
    
    fit <- fit_measles_largest
    cases_reported <- rstan::extract(fit, "v_cases_reported_sim")[[1]] %>% 
      as.data.frame() %>% 
      mutate(iteration=seq_along(V1)) %>% 
      pivot_longer(-iteration) %>% 
      group_by(iteration) %>% 
      mutate(time_reported=seq_along(iteration)-1) %>% 
      rename(cases_reported=value) %>% 
      ungroup() %>% 
      select(-name)
    
    cases_reported
  }),
  tar_target(plot_misfit_cases_per_day, {
    
    df_true <- measles_largest_cases %>% 
      mutate(cases_observed=c(0,diff(cases_reported)))
    
    measles_largest_sim_cases_observed %>% 
      ggplot() +
      geom_boxplot(aes(x=as.factor(time_reported), y=cases_reported)) +
      geom_point(data=df_true, aes(x=time_reported, y=cases_observed),
                colour="orange") +
      geom_line(data=df_true, aes(x=time_reported, y=cases_observed),
                colour="orange") +
      xlab("Reporting delay, days") +
      ylab("Cases reported on that day") +
      scale_x_discrete(breaks=seq(0, 40, 5)) +
      scale_y_continuous(breaks=seq(1, 10, 1))
  }),
  tar_target(plot_measles_misfit, {
    plot_grid(plot_measles_ecdf_absolute,
              plot_misfit_cases_per_day,
              labels=c("A.", "B."))
  }),
  tar_target(file_plot_measles_misfit, {
    filename <- "figures/measles_misfit.pdf"
    ggsave(filename,
           plot_measles_misfit,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_measles_misfit_png, {
    filename <- "figures/measles_misfit.png"
    ggsave(filename,
           plot_measles_misfit,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  
  # examining more dengue data in more depth
  tar_target(plot_dengue_ecdf_number, {
    
    df_summary <- dengue_processed %>% 
      filter(time_onset <= 156) %>% 
      group_by(time_onset) %>% 
      filter(max(cases_reported) > 5) %>% 
      mutate(
        ecdf=cases_reported/max(cases_reported),
        reporting_delay=seq_along(time_onset)
      )
    n_onset_periods <- 3
    df_sum <- df_summary %>% 
      ungroup() %>% 
      mutate(onset_period=case_when(
        time_onset<58~"A. early",
        time_onset>=58&time_onset<=108~"B. middle",
        time_onset>108~"C. late"
      )) %>% 
      group_by(time_onset) %>% 
      mutate(max_cases_reported=max(cases_reported)) %>% 
      ungroup()
    df_overall <- df_sum %>% 
      group_by(onset_period, reporting_delay) %>% 
      summarise(ecdf=mean(ecdf))
    onset_periods <- unique(df_overall$onset_period)
    d50s <- vector(length = n_onset_periods)
    for(i in 1:n_onset_periods) {
      df_tmp <- df_overall %>% filter(onset_period==onset_periods[i])
      f_iecdf <- approxfun(df_tmp$ecdf, df_tmp$reporting_delay)
      d50s[i] <- f_iecdf(0.5)
    }
    
    df_50s <- tibble(reporting_delay=d50s,
                     ecdf=0.5,
                     onset_period=onset_periods)
    
    df_sum %>% ggplot(aes(x=reporting_delay, y=ecdf)) +
      geom_line(alpha=0.4, aes(group=as.factor(time_onset),
                               colour=max_cases_reported)) +
      geom_line(data=df_overall, colour="orange", linetype=2) +
      geom_segment(data=df_50s, aes(xend=1, yend=ecdf),
                   colour="orange") +
      geom_segment(data=df_50s, aes(xend=reporting_delay, yend=0),
                   colour="orange") +
      facet_wrap(~onset_period) +
      scale_colour_viridis_c("Max\ncases", trans="sqrt") +
      scale_x_sqrt() +
      xlab("Reporting delay, weeks") +
      ylab("eCDF")
  }),
  tar_target(file_plot_dengue_ecdf_number, {
    filename <- "figures/dengue_ecdf.pdf"
    ggsave(filename,
           plot_dengue_ecdf_number,
           width = 10, height = 4);
    filename
  }),
  tar_target(file_plot_dengue_ecdf_number_png, {
    filename <- "figures/dengue_ecdf.png"
    ggsave(filename,
           plot_dengue_ecdf_number,
           width = 10, height = 4, dpi=1000);
    filename
  }),
  tar_target(dengue_delay_stan, {
    
    df_tmp <- dengue_processed %>% 
      group_by(time_onset) %>% 
      mutate(cases_true=max(cases_reported))
    
    reporting_days_sim <- seq(0, 40, 0.1)
    
    data_stan <- list(
      N_delay=nrow(df_tmp),
      time_reported=df_tmp$time_reported,
      time_onset=df_tmp$time_onset,
      cases_reported=df_tmp$cases_reported,
      cases_true=df_tmp$cases_true,
      n_reporting_window=1,
      reporting_window=rep(1, length(df_tmp$time_onset)),
      is_gamma_reporting_delay=0,
      reporting_days_sim=reporting_days_sim,
      N_sim=length(reporting_days_sim),
      is_binomial=1
    )
    
    options(mc.cores=4)
    init_fn <- function() {
      list(
        theta=matrix(c(10, 5), nrow = 1)
      )
    }
    optimizing(stan_model_delay, data=data_stan, init=init_fn, as_vector=FALSE)
  }),
  tar_target(dengue_delay_stan_ecdf, {
    fit <- dengue_delay_stan
    m_ecdf <- rstan::extract(fit, "m_ecdf_delay")[[1]][,,1]
    reporting_days_sim <- seq(0, 40, 0.1)
    lower <- apply(m_ecdf, 2, function(x) quantile(x, 0.025))
    upper <- apply(m_ecdf, 2, function(x) quantile(x, 0.975))
    middle <- apply(m_ecdf, 2, function(x) quantile(x, 0.5))
    df <- tibble(reporting_day=reporting_days_sim,
                 lower, middle, upper)
    df %>% 
      ggplot(aes(x=reporting_day, y=middle)) +
      geom_line() +
      geom_ribbon(aes(ymin=lower, ymax=upper))
    
    v_ecdf <- rstan::extract(fit, "v_ecdf")[[1]]
    
    df_tmp <- dengue_processed %>% 
      group_by(time_onset) %>% 
      mutate(cases_true=max(cases_reported))
    v_ecdf <- fit$par$v_ecdf
    v_cases_remaining <- fit$par$v_cases_remaining
    v_cases_observed <- fit$par$v_cases_observed
    v_p <- fit$par$v_p
    
    df_short <- tibble(
      time_onset=df_tmp$time_onset[1:(nrow(df_tmp) - 1)],
      time_reported=df_tmp$time_reported[1:(nrow(df_tmp) - 1)],
      cases_true=df_tmp$cases_true[1:(nrow(df_tmp) - 1)],
      cases_remaining=v_cases_remaining,
      cases_observed=v_cases_observed,
      p=v_p,
      ecdf=v_ecdf
    ) %>% 
      group_by(time_onset) %>% 
      mutate(reporting_delay=seq_along(time_onset)) %>% 
      filter(cases_remaining > 0) %>% 
      filter(cases_remaining != cases_observed) %>% # justified because each onset time always ends with one of these and we've filtered out the zero observations
      mutate(fraction_observed=cases_observed/cases_remaining)
    ggplot(df_short, aes(x=ecdf, y=fraction_observed)) +
      geom_jitter(aes(colour=cases_observed > 1), height = 0.02)
    hist(df_short$ecdf)
  }),
  
  # schematic figure
  tar_target(schematic_R, {
    c(rep(1.3, 40), rep(0.9, 20), seq(1, 1.2, length.out=40))
  }),
  tar_target(schematic_mu_delays, {
    seq(1, 20, length.out=20)
  }),
  tar_target(schematic_serial_parameters, {
    list(mean=5, sd=3)
  }),
  tar_target(schematic_cases, {
    days_total <- 100
    Rt_function <- stats::approxfun(1:days_total, schematic_R)
    s_params <- schematic_serial_parameters
    kappa <- 1000
    cases <- incidenceinflation::true_cases(days_total, Rt_function, kappa, s_params,
                                            initial_parameters=list(mean=30, length=30))
    df_true <- tibble(C=cases,
                      time_onset=seq_along(C))
    df_true
  }),
  tar_target(schematic_cases_reported_various_delays, {
    
    r_mus <- schematic_mu_delays
    r_params <- list(mean=3, sd=3)
    for(i in seq_along(r_mus)) {
      r_params$mean <- r_mus[i]
      df <- incidenceinflation::observed_cases(schematic_cases$C, r_params, days_max_follow_up=40) %>% 
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
    list(full=overall_df, snapshot=big_df)
  }),
  tar_target(schematic_model_poisson_rw_raw, "src/stan/poisson_rw.stan",
             format = "file"),
  tar_target(schematic_model_poisson_rw, stan_model(schematic_model_poisson_rw_raw)),
  tar_target(schematic_windows, unlist(map(seq(1, 100, 1), ~rep(., 1)))),
  tar_target(schematic_opt_R, {
    
    r_mus <- schematic_mu_delays
    s_params <- schematic_serial_parameters
    window <- schematic_windows
    overall_df <- schematic_cases_reported_various_delays$full
    
    for(i in seq_along(r_mus)) {

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
      
      opt <- rstan::optimizing(schematic_model_poisson_rw, data=data_stan, as_vector=FALSE)
      df_r_tmp <- tibble(R=opt$par$R,
                         index=seq(1, data_stan$K, 1),
                         delay=r_mus[i])
      
      if(i == 1)
        df_R <- df_r_tmp
      else
        df_R <- df_R %>% bind_rows(df_r_tmp)
    }
    df_R
  }),
  tar_target(schematic_lookup, {
    tibble(
      index=schematic_windows,
      time_onset=seq_along(index)
    )
  }),
  tar_target(schematic_t_threshold, 25),
  tar_target(schematic_R_cases, {
    
    df_R_processed <- schematic_opt_R %>% 
      left_join(schematic_lookup, by="index") %>% 
      filter(time_onset >= schematic_t_threshold) %>% 
      group_by(delay) %>% 
      mutate(time_onset=seq_along(R)) %>% 
      mutate(type="R[t]")
    # this is just for aesthetics
    t_before <- 46
    df_R_processed_1 <- df_R_processed %>% 
      filter(delay==max(df_R_processed$delay)) %>% 
      filter(time_onset < t_before)
    df_R_processed <- df_R_processed %>% 
      ungroup() %>% 
      filter(time_onset >= t_before) %>% 
      bind_rows(df_R_processed_1)
    
    df_cases_processed <- schematic_cases_reported_various_delays$snapshot %>% 
      filter(time_onset >= schematic_t_threshold) %>% 
      ungroup() %>% 
      group_by(delay) %>% 
      mutate(time_onset=seq_along(C)) %>% 
      mutate(type="Cases")
    
    df_cases_processed %>% 
      rename(value=C) %>% 
      bind_rows(df_R_processed %>% rename(value=R))
  }),
  tar_target(schematic_R_cases_true, {
    
    df_cases_true <- tibble(C=schematic_cases$C,
           time_onset=seq_along(C)) %>% 
      filter(time_onset >= schematic_t_threshold) %>% 
      mutate(time_onset=seq_along(C))
    
    df_R_true <- tibble(
      R=schematic_R,
      time_onset=seq_along(R)
    ) %>% 
      filter(time_onset >= schematic_t_threshold) %>% 
      mutate(time_onset=seq_along(R))
    
    df_cases_true %>% 
      mutate(type="Cases") %>% 
      rename(value=C) %>% 
      bind_rows(
        df_R_true %>% 
          mutate(type="R[t]") %>% 
          rename(value=R))
  }),
  tar_target(plot_schematic_top, {
    
    schematic_R_cases %>% 
      ggplot(aes(x=time_onset, y=value)) +
      geom_line(aes(colour=delay, group=delay)) +
      geom_line(data=schematic_R_cases_true, linetype=2) +
      scale_colour_viridis_c("Mean\ndelay", trans="reverse") +
      facet_wrap(~type, scales="free_y", labeller = label_parsed) +
      theme(
        axis.title.y = element_blank(),
        strip.text = element_text(size=18),
        text=element_text(size=14)
      ) +
      xlab("Onset day")
  }),
  tar_target(file_plot_schematic_top, {
    filename <- "figures/schematic_top.pdf"
    ggsave(filename,
           plot_schematic_top,
           width = 10, height = 4);
    filename
  }),
  tar_target(schematic_max_delay_df, {
    schematic_cases_reported_various_delays$full %>% 
      filter(delay==max(schematic_mu_delays))
  }),
  tar_target(schematic_moderate_delay_df, {
    schematic_cases_reported_various_delays$full %>% 
      filter(delay==5)
  }),
  tar_target(plot_schematic_trajectories, {
    
    df_tmp <- schematic_moderate_delay_df %>% 
      mutate(reporting_delay=time_reported-time_onset) %>% 
      filter(time_onset >= schematic_t_threshold) %>% 
      mutate(time_onset=time_onset - schematic_t_threshold + 1) %>% 
      mutate(time_reported=time_reported - schematic_t_threshold + 1)
    
    df_tmp %>%
      ggplot(aes(x=time_reported,
                 y=cases_reported,
                 group=as.factor(time_onset),
                 colour=reporting_delay)) +
      geom_line() +
      xlab("Onset day & report day") +
      ylab("Cases reported") +
      scale_colour_viridis_c("Delay", trans="reverse") +
      theme(legend.position = c(0.1, 0.74))
  }),
  tar_target(file_plot_schematic_trajectories, {
    filename <- "figures/schematic_trajectories.pdf"
    ggsave(filename,
           plot_schematic_trajectories,
           width = 6, height = 4);
    filename
  }),
  tar_target(schematic_rt_indices, {
    Rt_indices <- unlist(map(seq(1, 10, 1), ~rep(., 10)))
    lookup <- tibble(Rt_index=Rt_indices,
                     time_onset=seq_along(Rt_index))
    lookup
  }),
  tar_target(schematic_mcmc_fit, {
    
    # make problem easier by initialising near true values (this is only a schematic figure)
    inits_reporting <- list(mean=10, sd=3)
    lookup_true <- schematic_rt_indices %>% 
      mutate(Rt=schematic_R)
    inits_Rt <- lookup_true %>% 
      group_by(Rt_index) %>% 
      summarise(Rt=mean(Rt))
    inits_cases <- schematic_cases %>% 
      rename(cases_true=C)
    inits <- list(reporting=inits_reporting,
                  Rt=inits_Rt,
                  cases=inits_cases)
    fit_model_inits(schematic_moderate_delay_df %>% 
                select(-c(cases_true, delay)) %>% 
                left_join(schematic_rt_indices),
                schematic_serial_parameters,
                inits=inits,
                niterations=200)
    }),
  tar_target(plot_schematic_mcmc_cases, 
             plot_cases_simple_schematic(schematic_mcmc_fit$cases %>% 
                                   mutate(max_time="100"),
                                 schematic_moderate_delay_df,
                                 include_reporting=TRUE)),
  tar_target(file_plot_schematic_mcmc_cases, {
    filename <- "figures/schematic_cases.pdf"
    ggsave(filename,
           plot_schematic_mcmc_cases,
           width = 6, height = 4);
    filename
  }),
  tar_target(plot_schematic_rt, {
    
    summarise_schematic_rt <- schematic_mcmc_fit$Rt %>% 
      filter(iteration >= 100) %>% 
      group_by(Rt_index) %>% 
      summarise(
        lower=quantile(Rt, 0.025),
        middle=quantile(Rt, 0.5),
        upper=quantile(Rt, 0.975)
      )
    
    schematic_rt_df <- tibble(
      Rt=schematic_R,
      time_onset=seq_along(Rt)
    )
    
    schematic_rt_indices %>% 
      left_join(summarise_schematic_rt,
                by="Rt_index",
                relationship =
                  "many-to-many") %>%
      left_join(schematic_rt_df, by="time_onset") %>% 
      filter(time_onset>=25) %>% 
      mutate(time_onset=time_onset-25 +1) %>% 
      ggplot(aes(x=time_onset, y=middle)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4,
                  fill="#1B9E77") +
      geom_line(colour="#1B9E77") +
      geom_line(aes(y=Rt), linetype=2) +
      ylab(TeX("$R_t$")) +
      xlab("Onset day") +
      theme(legend.position="none") +
      scale_y_sqrt()
  }),
  tar_target(file_plot_schematic_rt, {
    filename <- "figures/schematic_rt.pdf"
    ggsave(filename,
           plot_schematic_rt,
           width = 6, height = 4);
    filename
  }),
  tar_target(reporting_parameters_5_3, {
    list(mean=5, sd=3)
  }),
  tar_target(plot_reporting_delay_schematic,
             plot_delay_true_series_simple(schematic_mcmc_fit$reporting %>% 
                                             mutate(max_time="100"),
                                           reporting_parameters_5_3)),
  tar_target(file_plot_reporting_delay_schematic, {
    filename <- "figures/schematic_reporting.pdf"
    ggsave(filename,
           plot_reporting_delay_schematic,
           width = 6, height = 4);
    filename
  }),
  
  # real Ebola data
  tar_target(ebola_data_real, {
    
    df <- read.csv("data/raw/evd_drc_2018-2020_daily.csv") %>% 
      rename(region=zds) %>% 
      mutate(date_onset=as.Date(date_onset, format="%Y-%m-%d")) %>% 
      filter(!is.na(date_onset)) %>% # there are 11 date_onsets that are NAs
      filter(!is.na(region)) # some NAs for region too
      
    date_min <- min(df$date_onset)
    date_max <- max(df$date_onset)
    
    regions <- unique(df$region)
    lookup <- expand_grid(
      date_onset=seq(date_min, date_max, 1),
      region=regions)
    df_c <- lookup %>% 
      left_join(df) %>% 
      mutate(n=if_else(is.na(n), 0, n))
    stopifnot(sum(df_c$n) == sum(df$n))
    
    df_c %>% 
      rename(cases_true=n)
  }),
  tar_target(ebola_data_beni_simulated, {
    # take last delay of 1.7 weeks = 11.9 days from https://www.sciencedirect.com/science/article/pii/S175543651830166X?via%3Dihub
    # assume a sd of 5 days
    
    r_params <- list(mean=11.9, sd=5)
    df_beni <- ebola_data_real %>% filter(region=="Beni")
    first_date <- df_beni$date_onset[which(df_beni$cases_true >= 1)[1]]
    df_beni <- df_beni %>%
      filter(date_onset >= first_date) %>% 
      mutate(time_onset=seq_along(date_onset))
    window_width <- 14
    n_rt <- floor(nrow(df_beni) / window_width)
    df_beni <- df_beni[1:(n_rt * window_width), ]
    Rt_indices <- unlist(map(seq(1, n_rt, 1), ~rep(., window_width)))
    df_beni <- df_beni %>% mutate(Rt_index=Rt_indices)
    lookup <- df_beni %>% select(date_onset, time_onset, Rt_index)
    df <- incidenceinflation::observed_cases(df_beni$cases_true, r_params, days_max_follow_up=40)
    
    df %>% 
      left_join(lookup)
  }),
  tar_target(ebola_data_beni_simulated_noisy, {
    # take last delay of 1.7 weeks = 11.9 days from https://www.sciencedirect.com/science/article/pii/S175543651830166X?via%3Dihub
    # assume a sd of 5 days
    
    r_params <- list(mean=11.9, sd=5)
    df_beni <- ebola_data_real %>% filter(region=="Beni")
    first_date <- df_beni$date_onset[which(df_beni$cases_true >= 1)[1]]
    df_beni <- df_beni %>%
      filter(date_onset >= first_date) %>% 
      mutate(time_onset=seq_along(date_onset))
    window_width <- 14
    n_rt <- floor(nrow(df_beni) / window_width)
    df_beni <- df_beni[1:(n_rt * window_width), ]
    Rt_indices <- unlist(map(seq(1, n_rt, 1), ~rep(., window_width)))
    df_beni <- df_beni %>% mutate(Rt_index=Rt_indices)
    lookup <- df_beni %>% select(date_onset, time_onset, Rt_index)
    df <- observed_cases(df_beni$cases_true, r_params, kappa=0.5, days_max_follow_up=40)
    
    df %>% 
      left_join(lookup)
  }),
  tar_target(fit_ebola_data_beni_simulated_310, {
    fit_ebola_real(310,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_330, {
    fit_ebola_real(330,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_350, {
    fit_ebola_real(350,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_370, {
    fit_ebola_real(370,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_390, {
    fit_ebola_real(390,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_410, {
    fit_ebola_real(410,
                   ebola_data_beni_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_310, {
    fit_ebola_real(310,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_330, {
    fit_ebola_real(330,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_350, {
    fit_ebola_real(350,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_370, {
    fit_ebola_real(370,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_390, {
    fit_ebola_real(390,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_beni_simulated_noisy_410, {
    fit_ebola_real(410,
                   ebola_data_beni_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_410, {
    
    time_end <- 410
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_390, {
    
    time_end <- 390
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_370, {
    
    time_end <- 370
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_350, {
    
    time_end <- 350
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_330, {
    
    time_end <- 330
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_310, {
    
    time_end <- 310
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_310, {
    
    time_end <- 310
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_330, {
    
    time_end <- 330
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_350, {
    
    time_end <- 350
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_370, {
    
    time_end <- 370
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_390, {
    
    time_end <- 390
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_beni_simulated_stan_noisy_410, {
    
    time_end <- 410
    
    fit_both_stan(
      model_poisson,
      ebola_data_beni_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_beni_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_beni_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(summarise_ebola_real_r_estimates, {
    
    # process stan fits
    R_1 <- summarise_r(fits_ebola_data_beni_simulated_stan_310) %>% 
      mutate(max_time=310)
    R_2 <- summarise_r(fits_ebola_data_beni_simulated_stan_330) %>% 
      mutate(max_time=330)
    R_3 <- summarise_r(fits_ebola_data_beni_simulated_stan_350) %>% 
      mutate(max_time=350)
    R_4 <- summarise_r(fits_ebola_data_beni_simulated_stan_370) %>% 
      mutate(max_time=370)
    R_5 <- summarise_r(fits_ebola_data_beni_simulated_stan_390) %>% 
      mutate(max_time=390)
    R_6 <- summarise_r(fits_ebola_data_beni_simulated_stan_410) %>% 
      mutate(max_time=410)
    R_stan <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_1 <- quantiles_r_other(fit_ebola_data_beni_simulated_310) %>% 
      mutate(max_time=310)
    R_2 <- quantiles_r_other(fit_ebola_data_beni_simulated_330) %>% 
      mutate(max_time=330)
    R_3 <- quantiles_r_other(fit_ebola_data_beni_simulated_350) %>% 
      mutate(max_time=350)
    R_4 <- quantiles_r_other(fit_ebola_data_beni_simulated_370) %>% 
      mutate(max_time=370)
    R_5 <- quantiles_r_other(fit_ebola_data_beni_simulated_390) %>% 
      mutate(max_time=390)
    R_6 <- quantiles_r_other(fit_ebola_data_beni_simulated_410) %>% 
      mutate(max_time=410)
    
    R_inflation <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(summarise_ebola_real_r_estimates_noisy, {
    
    # process stan fits
    R_1 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_310) %>% 
      mutate(max_time=310)
    R_2 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_330) %>% 
      mutate(max_time=330)
    R_3 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_350) %>% 
      mutate(max_time=350)
    R_4 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_370) %>% 
      mutate(max_time=370)
    R_5 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_390) %>% 
      mutate(max_time=390)
    R_6 <- summarise_r(fits_ebola_data_beni_simulated_stan_noisy_410) %>% 
      mutate(max_time=410)
    R_stan <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_1 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_310) %>% 
      mutate(max_time=310)
    R_2 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_330) %>% 
      mutate(max_time=330)
    R_3 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_350) %>% 
      mutate(max_time=350)
    R_4 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_370) %>% 
      mutate(max_time=370)
    R_5 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_390) %>% 
      mutate(max_time=390)
    R_6 <- quantiles_r_other(fit_ebola_data_beni_simulated_noisy_410) %>% 
      mutate(max_time=410)
    
    R_inflation <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(plot_ebola_real_rt, {
    
    lookup <- ebola_data_beni_simulated %>%
      select(time_onset, Rt_index, date_onset) %>% 
      unique()
    
    df <- lookup %>% 
      left_join(summarise_ebola_real_r_estimates) %>% 
    mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      filter(time_onset <= 410) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    
    df_obs <- df %>% 
      arrange(date_onset) %>% 
      group_by(max_time) %>% 
      summarise(observation_point=last(date_onset))
    
    df <- df %>% 
      left_join(df_obs)
    df %>% 
      filter(Rt_index >= 4) %>%  # initial period when highly uncertain
      filter(date_onset>="2019-02-01") %>% 
      filter(date_onset<="2019-09-10") %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~observation_point) +
      scale_y_sqrt() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
      geom_hline(yintercept = 1, linetype=2)
  }),
  tar_target(file_plot_ebola_real_rt, {
    filename <- "figures/ebola_real_rt.pdf"
    ggsave(filename,
           plot_ebola_real_rt,
           width = 10, height = 6);
    filename
  }),
  tar_target(plot_ebola_real_rt_noisy, {
    
    lookup <- ebola_data_beni_simulated %>%
      select(time_onset, Rt_index, date_onset) %>% 
      unique()
    
    df <- lookup %>% 
      left_join(summarise_ebola_real_r_estimates_noisy) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      filter(time_onset <= 410) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    
    df_obs <- df %>% 
      arrange(date_onset) %>% 
      group_by(max_time) %>% 
      summarise(observation_point=last(date_onset))
    
    df <- df %>% 
      left_join(df_obs)
    df %>% 
      filter(Rt_index >= 4) %>%  # initial period when highly uncertain
      filter(date_onset>="2019-02-01") %>% 
      filter(date_onset<="2019-09-10") %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~observation_point) +
      scale_y_sqrt() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
      geom_hline(yintercept = 1, linetype=2)
  }),
  tar_target(file_plot_ebola_real_rt_noisy, {
    filename <- "figures/ebola_real_rt_noisy.pdf"
    ggsave(filename,
           plot_ebola_real_rt_noisy,
           width = 10, height = 6);
    filename
  }),
  tar_target(ebola_data_all_simulated, {
    # take last delay of 1.7 weeks = 11.9 days from https://www.sciencedirect.com/science/article/pii/S175543651830166X?via%3Dihub
    # assume a sd of 5 days
    
    r_params <- list(mean=11.9, sd=5)
    df_beni <- ebola_data_real %>% 
      group_by(date_onset) %>% 
      summarise(cases_true=sum(cases_true)) %>% 
      ungroup()
    first_date <- df_beni$date_onset[which(df_beni$cases_true >= 1)[1]]
    df_beni <- df_beni %>%
      filter(date_onset >= first_date) %>% 
      mutate(time_onset=seq_along(date_onset))
    window_width <- 14
    n_rt <- floor(nrow(df_beni) / window_width)
    df_beni <- df_beni[1:(n_rt * window_width), ]
    Rt_indices <- unlist(map(seq(1, n_rt, 1), ~rep(., window_width)))
    df_beni <- df_beni %>% mutate(Rt_index=Rt_indices)
    lookup <- df_beni %>% select(date_onset, time_onset, Rt_index)
    df <- incidenceinflation::observed_cases(df_beni$cases_true, r_params, days_max_follow_up=40)
    
    df %>% 
      left_join(lookup)
  }),
  tar_target(ebola_data_all_simulated_noisy, {
    # take last delay of 1.7 weeks = 11.9 days from https://www.sciencedirect.com/science/article/pii/S175543651830166X?via%3Dihub
    # assume a sd of 5 days
    
    r_params <- list(mean=11.9, sd=5)
    df_beni <- ebola_data_real %>% 
      group_by(date_onset) %>% 
      summarise(cases_true=sum(cases_true)) %>% 
      ungroup()
    first_date <- df_beni$date_onset[which(df_beni$cases_true >= 1)[1]]
    df_beni <- df_beni %>%
      filter(date_onset >= first_date) %>% 
      mutate(time_onset=seq_along(date_onset))
    window_width <- 14
    n_rt <- floor(nrow(df_beni) / window_width)
    df_beni <- df_beni[1:(n_rt * window_width), ]
    Rt_indices <- unlist(map(seq(1, n_rt, 1), ~rep(., window_width)))
    df_beni <- df_beni %>% mutate(Rt_index=Rt_indices)
    lookup <- df_beni %>% select(date_onset, time_onset, Rt_index)
    df <- observed_cases(df_beni$cases_true, r_params, kappa=0.5, days_max_follow_up=40)
    
    df %>% 
      left_join(lookup)
  }),
  tar_target(fit_ebola_data_all_simulated_310, {
    fit_ebola_real(310,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_330, {
    fit_ebola_real(330,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_350, {
    fit_ebola_real(350,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }), 
  tar_target(fit_ebola_data_all_simulated_370, {
    fit_ebola_real(370,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_390, {
    fit_ebola_real(390,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_410, {
    fit_ebola_real(410,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_430, {
    fit_ebola_real(430,
                   ebola_data_all_simulated,
                   serial_parameters_ebola,
                   niterations=50)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_310, {
    fit_ebola_real(310,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_330, {
    fit_ebola_real(330,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_350, {
    fit_ebola_real(350,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_370, {
    fit_ebola_real(370,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_390, {
    fit_ebola_real(390,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_410, {
    fit_ebola_real(410,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=400)
  }),
  tar_target(fit_ebola_data_all_simulated_noisy_430, {
    fit_ebola_real(430,
                   ebola_data_all_simulated_noisy,
                   serial_parameters_ebola,
                   niterations=50)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_430, {
    
    time_end <- 430
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_410, {
    
    time_end <- 410
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_390, {
    
    time_end <- 390
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_370, {
    
    time_end <- 370
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_350, {
    
    time_end <- 350
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_330, {
    
    time_end <- 330
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_310, {
    
    time_end <- 310
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_310, {
    
    time_end <- 310
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_330, {
    
    time_end <- 330
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_350, {
    
    time_end <- 350
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_370, {
    
    time_end <- 370
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_390, {
    
    time_end <- 390
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_410, {
    
    time_end <- 410
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(fits_ebola_data_all_simulated_stan_noisy_430, {
    
    time_end <- 430
    
    fit_both_stan(
      model_poisson,
      ebola_data_all_simulated_noisy %>%
        select(-c(cases_true, date_onset)) %>% 
        filter(time_onset <= time_end) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      ebola_data_all_simulated_noisy %>% 
        filter(time_reported <= time_end) %>% 
        select(-c(cases_true, date_onset)) %>% 
        group_by(time_onset) %>% 
        summarise(C=last(cases_reported)),
      40,
      ebola_generation_ws,
      ebola_data_all_simulated_noisy %>% 
        filter(time_onset <= time_end) %>% 
        select(time_onset, Rt_index) %>% 
        unique() %>% 
        pull(Rt_index),
      I_init=rep(1, 40),
      n_iterations=2000)
  }),
  tar_target(summarise_ebola_real_r_estimates_all, {
    
    # process stan fits
    R_1 <- summarise_r(fits_ebola_data_all_simulated_stan_310) %>% 
      mutate(max_time=310)
    R_2 <- summarise_r(fits_ebola_data_all_simulated_stan_330) %>% 
      mutate(max_time=330)
    R_3 <- summarise_r(fits_ebola_data_all_simulated_stan_350) %>% 
      mutate(max_time=350)
    R_4 <- summarise_r(fits_ebola_data_all_simulated_stan_370) %>% 
      mutate(max_time=370)
    R_5 <- summarise_r(fits_ebola_data_all_simulated_stan_390) %>% 
      mutate(max_time=390)
    R_6 <- summarise_r(fits_ebola_data_all_simulated_stan_410) %>% 
      mutate(max_time=410)
    R_7 <- summarise_r(fits_ebola_data_all_simulated_stan_410) %>% 
      mutate(max_time=430)
    R_stan <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6,
        R_7
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_1 <- quantiles_r_other(fit_ebola_data_all_simulated_310) %>% 
      mutate(max_time=310)
    R_2 <- quantiles_r_other(fit_ebola_data_all_simulated_330) %>% 
      mutate(max_time=330)
    R_3 <- quantiles_r_other(fit_ebola_data_all_simulated_350) %>% 
      mutate(max_time=350)
    R_4 <- quantiles_r_other(fit_ebola_data_all_simulated_370) %>% 
      mutate(max_time=370)
    R_5 <- quantiles_r_other(fit_ebola_data_all_simulated_390) %>% 
      mutate(max_time=390)
    R_6 <- quantiles_r_other(fit_ebola_data_all_simulated_410) %>% 
      mutate(max_time=410)
    R_7 <- quantiles_r_other(fit_ebola_data_all_simulated_430) %>% 
      mutate(max_time=430)
    
    R_inflation <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6,
        R_7
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(summarise_ebola_real_r_estimates_noisy_all, {
    
    # process stan fits
    R_1 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_310) %>% 
      mutate(max_time=310)
    R_2 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_330) %>% 
      mutate(max_time=330)
    R_3 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_350) %>% 
      mutate(max_time=350)
    R_4 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_370) %>% 
      mutate(max_time=370)
    R_5 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_390) %>% 
      mutate(max_time=390)
    R_6 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_410) %>% 
      mutate(max_time=410)
    R_7 <- summarise_r(fits_ebola_data_all_simulated_stan_noisy_430) %>% 
      mutate(max_time=430)
    R_stan <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6,
        R_7
      ) %>% 
      mutate(estimator="stan")
    
    # process other fits
    R_1 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_310) %>% 
      mutate(max_time=310)
    R_2 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_330) %>% 
      mutate(max_time=330)
    R_3 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_350) %>% 
      mutate(max_time=350)
    R_4 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_370) %>% 
      mutate(max_time=370)
    R_5 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_390) %>% 
      mutate(max_time=390)
    R_6 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_410) %>% 
      mutate(max_time=410)
    R_7 <- quantiles_r_other(fit_ebola_data_all_simulated_noisy_430) %>% 
      mutate(max_time=430)
    
    R_inflation <- R_1 %>% 
      bind_rows(
        R_2,
        R_3,
        R_4,
        R_5,
        R_6,
        R_7
      ) %>% 
      mutate(estimator="inflation")
    R_stan %>% 
      bind_rows(R_inflation)
  }),
  tar_target(data_plot_ebola_real_rt_all, {
    lookup <- ebola_data_all_simulated %>%
      select(time_onset, Rt_index, date_onset) %>% 
      unique()
    
    df <- lookup %>% 
      left_join(summarise_ebola_real_r_estimates_all) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      filter(time_onset <= 406) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    df_obs <- df %>% 
      arrange(date_onset) %>% 
      group_by(max_time) %>% 
      summarise(observation_point=last(date_onset))
    
    df <- df %>% 
      left_join(df_obs) %>% 
      mutate(observation_point=format(observation_point, "%d %b %Y"))
    df %>% 
      filter(Rt_index >= 4) %>%  # initial period when highly uncertain
      filter(date_onset>="2018-10-01") %>% 
      filter(max_time %in% c(330, 350, 370, 390, 410))
  }),
  tar_target(plot_ebola_real_rt_all, {
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    
    data_plot_ebola_real_rt_all %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~observation_point) +
      scale_y_continuous(limits=c(0, 4)) +
      scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
      geom_hline(yintercept = 1, linetype=2) +
      theme(legend.position = c(0.8, 0.25))
  }),
  tar_target(file_plot_ebola_real_rt_all, {
    filename <- "figures/ebola_real_rt_all.pdf"
    ggsave(filename,
           plot_ebola_real_rt_all,
           width = 10, height = 6);
    filename
  }),
  tar_target(data_plot_ebola_real_rt_noisy_all, {
    
    lookup <- ebola_data_all_simulated %>%
      select(time_onset, Rt_index, date_onset) %>% 
      unique()
    
    df <- lookup %>% 
      left_join(summarise_ebola_real_r_estimates_noisy_all) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      filter(time_onset <= 406) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    df_obs <- df %>% 
      arrange(date_onset) %>% 
      group_by(max_time) %>% 
      summarise(observation_point=last(date_onset))
    df <- df %>% 
      left_join(df_obs) %>% 
      mutate(observation_point=format(observation_point, "%d %b %Y"))
    df %>% 
      filter(Rt_index >= 4) %>%  # initial period when highly uncertain
      filter(date_onset>="2018-10-01") %>% 
      filter(max_time %in% c(330, 350, 370, 390, 410))
  }),
  tar_target(plot_ebola_real_rt_noisy_all, {
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    
    data_plot_ebola_real_rt_noisy_all %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_wrap(~observation_point) +
      scale_y_continuous(limits=c(0, 4)) +
      scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
      geom_hline(yintercept = 1, linetype=2) +
      theme(legend.position = c(0.8, 0.25))
  }),
  tar_target(file_plot_ebola_real_rt_noisy_all, {
    filename <- "figures/ebola_real_rt_noisy_all.pdf"
    ggsave(filename,
           plot_ebola_real_rt_noisy_all,
           width = 10, height = 6);
    filename
  }),
  tar_target(plot_ebola_real_rt_both, {
    
    cols <- c(naive="#FF8C00", "gold standard"="black", incidenceinflation="#1B9E77")
    
    data_plot_ebola_real_rt_all %>% 
      mutate(type="A. correct") %>% 
      bind_rows(data_plot_ebola_real_rt_noisy_all %>% 
                  mutate(type="B. misspecified")) %>% 
      ggplot(aes(x=date_onset, y=middle, fill=variant)) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
      geom_line(aes(colour=variant)) +
      scale_color_manual("Estimator", values=cols) +
      scale_fill_manual("Estimator", values=cols) +
      ylab(TeX("$R_t$")) +
      xlab("Onset date") +
      facet_grid(vars(observation_point), vars(type), scales="free") +
      scale_y_continuous() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
      geom_hline(yintercept = 1, linetype=2) +
      theme(text=element_text(size=16))
  }),
  tar_target(file_plot_ebola_real_rt_both, {
    filename <- "figures/ebola_real_rt_both.pdf"
    ggsave(filename,
           plot_ebola_real_rt_both,
           width = 10, height = 8);
    filename
  }),
  
  tar_target(ebola_real_combined, {
    
    df_1 <- fit_ebola_data_all_simulated_410$cases %>% 
      mutate(max_time="410")
    df_2 <- fit_ebola_data_all_simulated_390$cases %>% 
      mutate(max_time="390")
    df_3 <- fit_ebola_data_all_simulated_370$cases %>% 
      mutate(max_time="370")
    df_4 <- fit_ebola_data_all_simulated_350$cases %>% 
      mutate(max_time="350")
    df_5 <- fit_ebola_data_all_simulated_330$cases %>% 
      mutate(max_time="330")
    df_6 <- fit_ebola_data_all_simulated_310$cases %>% 
      mutate(max_time="310")
    
    df_1 %>% 
      bind_rows(
        df_2,
        df_3,
        df_4,
        df_5
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "410", "390", "370", "350", "330")) %>% 
      mutate(max_time=fct_rev(max_time)) %>% 
      filter(time_onset <= 406) # to match Rt
  }),
  tar_target(ebola_real_noisy_combined, {
    
    df_1 <- fit_ebola_data_all_simulated_noisy_410$cases %>% 
      mutate(max_time="410")
    df_2 <- fit_ebola_data_all_simulated_noisy_390$cases %>% 
      mutate(max_time="390")
    df_3 <- fit_ebola_data_all_simulated_noisy_370$cases %>% 
      mutate(max_time="370")
    df_4 <- fit_ebola_data_all_simulated_noisy_350$cases %>% 
      mutate(max_time="350")
    df_5 <- fit_ebola_data_all_simulated_noisy_330$cases %>% 
      mutate(max_time="330")
    df_6 <- fit_ebola_data_all_simulated_noisy_310$cases %>% 
      mutate(max_time="310")
    
    df_1 %>% 
      bind_rows(
        df_2,
        df_3,
        df_4,
        df_5
      ) %>% 
      mutate(max_time=as.factor(as.character(max_time))) %>% 
      mutate(max_time=fct_relevel(max_time, "410", "390", "370", "350", "330")) %>% 
      mutate(max_time=fct_rev(max_time)) %>% 
      filter(time_onset <= 406) # to match Rt
  }),
  tar_target(cases_ebola_real_true, {
    
    d <- ebola_data_all_simulated %>% 
      group_by(date_onset) %>% 
      summarise(cases_true=first(cases_true))
    
    e <- tribble(
      ~max_time,
      410,
      390,
      370,
      350,
      330)
    
    for(i in seq_along(e$max_time)) {
      tmp <- ebola_data_all_simulated %>% 
        mutate(max_time=e$max_time[i]) %>% 
        filter(time_reported <= max_time)
      if(i == 1)
        big_df <- tmp
      else
        big_df <- big_df %>% bind_rows(tmp)
    }
    df_reported <- big_df %>% 
      group_by(max_time, time_onset, date_onset) %>% 
      summarise(cases_reported=last(cases_reported))
    df_reported %>%
      left_join(d, by="date_onset") %>% 
      mutate(max_time=as.factor(as.character(max_time)))
  }),
  tar_target(ebola_observation_dates, {
    
    lookup <- ebola_data_all_simulated %>%
      select(time_onset, Rt_index, date_onset) %>% 
      unique()
    
    df <- lookup %>% 
      left_join(summarise_ebola_real_r_estimates_noisy_all) %>% 
      mutate(variant=paste0(estimator, "-", type)) %>% 
      mutate(max_time=as.factor(max_time)) %>% 
      filter(time_onset <= 406) %>% 
      mutate(variant=case_when(
        variant=="stan-full"~"gold standard",
        variant=="stan-snapshot"~"naive",
        variant=="inflation-snapshot"~"incidenceinflation"
      )) %>% 
      mutate(variant=as.factor(variant)) %>% 
      mutate(variant=fct_relevel(variant, "naive", "gold standard", "incidenceinflation"))
    
    df %>% 
      arrange(date_onset) %>% 
      group_by(max_time) %>% 
      summarise(observation_point=last(date_onset))
  }),
  tar_target(plot_cases_ebola_combined,
             plot_cases_vs_true_series_3(ebola_real_combined,
                                         cases_ebola_real_true,
                                         ebola_data_all_simulated %>%
                                           filter(date_onset>="2019-01-01"),
                                         observation_dates=ebola_observation_dates,
                                         transform_type="sqrt",
                                         include_reporting=TRUE,
                                         ylimits=c(0, 50),
                                         min_date="2019-01-01")),
  tar_target(file_plot_cases_ebola_combined, {
    filename <- "figures/ebola_real_cases.pdf"
    ggsave(filename,
           plot_cases_ebola_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(plot_cases_ebola_noisy_combined,
             plot_cases_vs_true_series_3(ebola_real_noisy_combined,
                                         cases_ebola_real_true,
                                         ebola_data_all_simulated %>%
                                           filter(date_onset>="2019-01-01"),
                                         observation_dates=ebola_observation_dates,
                                         transform_type="sqrt",
                                         include_reporting=TRUE,
                                         ylimits=c(0, 50),
                                         min_date="2019-01-01")),
  tar_target(file_plot_cases_ebola_noisy_combined, {
    filename <- "figures/ebola_real_noisy_cases.pdf"
    ggsave(filename,
           plot_cases_ebola_noisy_combined,
           width = 10, height = 4);
    filename
  }),
  tar_target(plot_cases_ebola_real_both, {
    df_correct <- process_cases_data(ebola_real_combined,
                                     cases_ebola_real_true,
                                     ebola_data_all_simulated %>%
                                       filter(date_onset>="2019-01-01"),
                                     observation_dates=ebola_observation_dates) %>% 
      mutate(type="A. correct")
    df_mis <- process_cases_data(ebola_real_noisy_combined,
                                     cases_ebola_real_true,
                                     ebola_data_all_simulated %>%
                                       filter(date_onset>="2019-01-01"),
                                     observation_dates=ebola_observation_dates) %>% 
      mutate(type="B. misspecified")
    
    df_both <- df_correct %>% 
      bind_rows(df_mis)
    
    min_date <- "2019-01-01"
    df_both <- df_both %>% filter(date_onset>=min_date)
      
    df_both %>% 
        ggplot(aes(x=date_onset)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="#1B9E77", alpha=0.4) +
        geom_line(aes(y=cases_reported), colour="#FF8C00") +
        geom_line(aes(y=middle), colour="#1B9E77") +
        geom_line(aes(y=cases_true), colour="black") +
        xlab("Onset date") +
        ylab("Cases") +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
        facet_grid(vars(observation_date), vars(type), scales="free") +
        scale_y_sqrt() +
        theme(text=element_text(size=16))
    }),
  tar_target(file_plot_cases_ebola_real_both, {
    filename <- "figures/ebola_real_cases_both.pdf"
    ggsave(filename,
           plot_cases_ebola_real_both,
           width = 10, height = 8);
    filename
  }),
  tar_target(plot_ebola_real_incidence, {
    
    ebola_data_real %>%
      group_by(date_onset) %>%
      summarise(cases_true=sum(cases_true)) %>%
      ggplot(aes(x=date_onset, y=cases_true)) +
      geom_line()
  })
)
