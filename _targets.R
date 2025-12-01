library(targets)
library(crew)

# Set up controller to launch 3 workers on local processes
controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 3,
  seconds_idle = 10,
  crashes_error = 10
)

# Set targets options and required packages
tar_option_set(packages = c("tibble", "readr", "dplyr", "tidyr", "purrr", "rlang", 
                            "ggplot2", "ggdist", "patchwork", "ggtext",
                            "tRophicPosition", "MixSIAR", "compositions",
                            "brms", "cmdstanr", "tidybayes"),
               controller = controller,
               memory = "transient", 
               garbage_collection = TRUE,
               seed = 25)

# Load functions
targets::tar_source()

# Set number of simulations
nrep <- 100

# Pipeline
list(
  # Simulations with one or two baselines 
  tar_target(sim_baselines, run_simulations_baselines(nrep = nrep, 
                                                      C_baselines = 20,
                                                      C_consumers = 40,
                                                      run = "short")),
  tar_target(performance_pk, get_performance_pk(sim_baselines, nrep = nrep)),
  tar_target(performance_TP, get_performance_TP(sim_baselines, nrep = nrep)),
  tar_target(model_baseline, fit_baseline_effect(performance_pk, performance_TP)),
  tar_target(baseline_contrasts, get_baseline_contrasts(model_baseline)),
  tar_target(plot_baseline, plot_baseline_effect(performance_pk,
                                                 performance_TP,
                                                 model_baseline,
                                                 baseline_contrasts)),
  tar_target(plot_pk_vs_tp, plot_TP_error_effect(sim_baselines,
                                                 performance_pk,
                                                 performance_TP,
                                                 model_baseline,
                                                 .width = ppoints(50),
                                                 alpha_points = 0.2)),
  
  # Simulations with different consumer sample sizes
  tar_target(sim_food_web, simulate_food_web(n_sources = 3,
                                             n_baselines = 2,
                                             n_families = 1,
                                             n_species = 1,
                                             C_baselines = 20,
                                             C_consumers = 40)),
  tar_target(sim_sample_size_consumer, run_simulations_consumer_sample_size(sim_food_web,
                                                                            sample_size = c(1, 3, 5, 10, 20), 
                                                                            nrep = nrep, 
                                                                            C = 40,
                                                                            run = "short")),
  tar_target(model_median_sample_size_consumer, fit_sample_size_effect(sim_sample_size_consumer,
                                                                       response = "alpha_median")),
  tar_target(slopes_median_sample_size_consumer, get_slopes(model_median_sample_size_consumer)),
  tar_target(model_CI_sample_size_consumer, fit_sample_size_effect(sim_sample_size_consumer, 
                                                                   response = "alpha_CI")),
  tar_target(slopes_CI_sample_size_consumer, get_slopes(model_CI_sample_size_consumer)),
  tar_target(plot_sample_size_consumer, plot_sample_size_effect(sim_food_web, 
                                                                sim_sample_size_consumer, 
                                                                model_median_sample_size_consumer, 
                                                                slopes_median_sample_size_consumer,
                                                                model_CI_sample_size_consumer,
                                                                slopes_CI_sample_size_consumer,
                                                                xlab = "Consumer sample size",
                                                                .width = ppoints(50),
                                                                alpha_points = 0.2)),
  
  # Simulations with different sources sample sizes
  tar_target(sim_sample_size_source, run_simulations_source_sample_size(sim_food_web,
                                                                        sample_size = c(1, 3, 5, 10, 20), 
                                                                        nrep = nrep,
                                                                        run = "short")),
  tar_target(model_median_sample_size_source, fit_sample_size_effect(sim_sample_size_source,
                                                                     response = "alpha_median")),
  tar_target(slopes_median_sample_size_source, get_slopes(model_median_sample_size_source)),
  tar_target(model_CI_sample_size_source, fit_sample_size_effect(sim_sample_size_source, 
                                                                 response = "alpha_CI")),
  tar_target(slopes_CI_sample_size_source, get_slopes(model_CI_sample_size_source)),
  tar_target(plot_sample_size_source, plot_sample_size_effect(sim_food_web, 
                                                              sim_sample_size_source, 
                                                              model_median_sample_size_source, 
                                                              slopes_median_sample_size_source,
                                                              model_CI_sample_size_source,
                                                              slopes_CI_sample_size_source,
                                                              xlab = "Source sample size",
                                                              .width = ppoints(50),
                                                              alpha_points = 0.2))
)
