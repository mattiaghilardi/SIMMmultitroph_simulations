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
tar_source()

# Number of replicates for the baseline and sample-size simulations
# (the TDF uncertainty simulations use nrep_tdf() below)
nrep <- 100

# TDF uncertainty simulations ------------------------------------------------
# One target per simulation and consumer sample size, so that each is saved
# when it finishes and targets can run in parallel on the crew workers.
# Data are seeded per replicate inside the functions, so results do not depend
# on how the simulations are split into targets.

# Number of replicates by consumer sample size (more at n = 1, where estimates vary most)
nrep_tdf <- function(n) if (n == 1) 50 else 25
sample_sizes_tdf <- c(1, 5, 10, 20)
tdf_sd_levels_fit <- c("zero", "one_step", "scaled_TP")

# Error structures by sample size (process only at n = 1)
err_for <- function(n) {
  if (n == 1) "process"
  else if (n == 20) "process_residual"
  else c("process_residual", "process")
}

# Combine targets with dplyr::bind_rows()
combine_target <- function(name, parts) {
  tar_target_raw(name, as.call(c(quote(dplyr::bind_rows), purrr::map(parts, as.symbol))))
}

# TDF SD supplied to MixSIAR, without or with individual variation in TDF
grid_tdf_sd <- expand.grid(tdf_variation = c("none", "per_individual"),
                           n = sample_sizes_tdf,
                           stringsAsFactors = FALSE)

targets_tdf_sd <- purrr::map2(grid_tdf_sd$tdf_variation, grid_tdf_sd$n, function(v, n) {
  tar_target_raw(
    name = paste("sim_tdf_sd", v, n, sep = "_"),
    command = bquote(run_simulations_tdf_sd(sim_food_web,
                                            sample_size = .(n),
                                            nrep = .(nrep_tdf(n)),
                                            C = 40,
                                            tdf_variation = .(v),
                                            error_structures = .(err_for(n)),
                                            tdf_sd_levels = .(tdf_sd_levels_fit),
                                            run = "short"))
  )
})

# Misspecified mean TDF
targets_tdf_offset <- purrr::map(sample_sizes_tdf, function(n) {
  tar_target_raw(
    name = paste("sim_tdf_offset", n, sep = "_"),
    command = bquote(run_simulations_tdf_offset(sim_food_web,
                                                sample_size = .(n),
                                                nrep = .(nrep_tdf(n)),
                                                C = 40,
                                                tdf_variation = "none",
                                                tdf_sd_levels = .(tdf_sd_levels_fit),
                                                error_structures = "process_residual",
                                                run = "short"))
  )
})

targets_tdf_combined <- list(
  combine_target("sim_tdf_sd_none", paste0("sim_tdf_sd_none_", sample_sizes_tdf)),
  combine_target("sim_tdf_sd_per_individual", paste0("sim_tdf_sd_per_individual_", sample_sizes_tdf)),
  combine_target("sim_tdf_offset_all", paste0("sim_tdf_offset_", sample_sizes_tdf))
)

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
  tar_target(coverage_tp_error, get_coverage_tp_error(performance_pk, 
                                                      performance_TP, 
                                                      n_bins = 5)),
  tar_target(model_coverage_tp_error, fit_coverage_tp_error(performance_pk, performance_TP)),
  tar_target(plot_pk_vs_tp, plot_TP_error_effect(performance_pk,
                                                 performance_TP,
                                                 model_baseline,
                                                 coverage_tp_error,
                                                 model_coverage_tp_error,
                                                 .width = ppoints(50),
                                                 alpha_points = 0.2)),
  tar_target(plot_tp_error_analytic, plot_TP_error_analytic(sim_baselines,
                                                            performance_TP,
                                                            alpha_points = 0.2)),
  tar_target(baseline_contrasts,
             get_baseline_contrasts(model_baseline, model_coverage_tp_error)),
  tar_target(baseline_slopes,
             list(max_error = get_baseline_slopes(model_baseline),
                  coverage = get_baseline_slopes(model_coverage_tp_error))),
  tar_target(baseline_contrast_curves,
             get_baseline_contrast_curves(performance_TP, model_baseline,
                                          model_coverage_tp_error)),
  tar_target(plot_contrasts_baselines,
             plot_baseline_contrasts(baseline_contrast_curves, 
                                     baseline_contrasts,
                                     .width = ppoints(50))),
  tar_target(coverage_contrast_summary,
             summarise_coverage_contrast(baseline_contrast_curves)),
  
  # Simulations with different consumer sample sizes
  tar_target(sim_food_web, simulate_food_web(n_sources = 3,
                                             n_baselines = 2,
                                             n_families = 1,
                                             n_species = 1,
                                             C_baselines = 20,
                                             C_consumers = 40)),
  tar_target(limits_sample_size, get_sample_size_limits(sim_food_web)),
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
                                                                alpha_points = 0.2,
                                                                limits = dplyr::filter(limits_sample_size, 
                                                                                       simulation == "consumer"))),
  
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
                                                              alpha_points = 0.2,
                                                              limits = dplyr::filter(limits_sample_size, 
                                                                                     simulation == "source"))),
  
  # Simulations with TDF uncertainty
  targets_tdf_sd,
  targets_tdf_offset,
  targets_tdf_combined,
  tar_target(performance_tdf_sd, summarise_tdf_sd(dplyr::bind_rows(sim_tdf_sd_none,
                                                                   sim_tdf_sd_per_individual))),
  tar_target(plot_tdf_sd_none, plot_tdf_sd(sim_tdf_sd_none,
                                           limits = dplyr::filter(limits_sample_size,
                                                                  simulation == "consumer"))),
  tar_target(plot_tdf_sd_per_individual, plot_tdf_sd(sim_tdf_sd_per_individual,
                                                     limits = dplyr::filter(limits_sample_size,
                                                                            simulation == "consumer"))),
  tar_target(performance_tdf_offset, summarise_tdf_offset(sim_tdf_offset_all)),
  tar_target(plot_sim_tdf_offset, plot_tdf_offset(sim_tdf_offset_all, sim_food_web))
)
