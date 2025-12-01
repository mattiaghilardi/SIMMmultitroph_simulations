#' Test effect of consumers sample size on mixing model accuracy
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param sample_size A numeric vector of sample sizes to simulate
#' @param nrep Number of simulations to run
#' @param C Total consumer consumption rate. Default to 20
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A list
run_simulations_consumer_sample_size <- function(sim_food_web,
                                                 sample_size = c(1, 3, 5, 10, 20), 
                                                 nrep = 100,
                                                 C = 20,
                                                 run = "test") {
  
  # Simulate "nrep" consumer data sets for a range of sample sizes
  # using fixed sources, baselines, alpha and TP of consumer
  # Estimate source contribution for each and return summary table
  sim_sample_size <- purrr::map(
    sample_size |> 
      rlang::set_names(),
    function(i) {
      purrr::map(
        1:nrep,
        function(j) {
          # Simulate data
          data <- simulate_consumer_data(i, 
                                         sources = sim_food_web$sources_summary,
                                         alpha = sim_food_web$alpha_species$alpha,
                                         C = C,
                                         TP_mean = sim_food_web$TP_species$TP, 
                                         TP_sd = 0.05,
                                         deltaC = sim_food_web$TDF$deltaC, 
                                         deltaN = sim_food_web$TDF$deltaN, 
                                         seed = j)
          
          # Estimate TP
          si_data_list <- tRophicPosition::loadIsotopeData(sim_food_web$isotopes |> 
                                                             filter(type == "baseline") |> 
                                                             bind_rows(data |> 
                                                                         mutate(type = "higher consumer",
                                                                                name = "species 1")),
                                                           consumer = "species 1",
                                                           b1 = "primary consumer 1", 
                                                           b2 = "primary consumer 2", 
                                                           baselineColumn = "name", 
                                                           consumersColumn = "name",
                                                           groupsColumn = NULL,
                                                           d13C = "d13C",
                                                           d15N = "d15N",
                                                           deltaC = sim_food_web$TDF$deltaC,
                                                           deltaN = sim_food_web$TDF$deltaN)
          
          model.string <- tRophicPosition::jagsBayesianModel(model = "twoBaselinesFull", lambda = 2)
          
          fit_TP <- tRophicPosition::TPmodel(data = si_data_list, 
                                             model.string = model.string, 
                                             n.adapt = 2000, 
                                             n.chains = 4, 
                                             quiet = FALSE)
          
          post <- tRophicPosition::posteriorTP(fit_TP,
                                               variable.names = c("TP", "alpha"),
                                               n.iter = 2000,
                                               burnin = 2000,
                                               thin = 1,
                                               quiet = FALSE)
          
          post_combined <- bind_rows(purrr::map(post, as.data.frame))
          
          TP <- tRophicPosition::getPosteriorMode(post_combined)$TP
          
          # Rescale consumers
          mean_TDF_C <- mean(sim_food_web$TDF$deltaC)
          mean_TDF_N <- mean(sim_food_web$TDF$deltaN)
          
          consumers_rescaled <- data |> 
            mutate(d13C = d13C - mean_TDF_C * (TP - 1),
                   d15N = d15N - mean_TDF_N * (TP - 1))
          
          # Mixing model
          sourcespath <- tempfile(pattern = paste("sim_sources", i, j, sep = "_"), fileext = ".csv")
          write_csv(sim_food_web$isotopes |> 
                      filter(type == "source") |> 
                      select(source = name, d13C, d15N), 
                    file = sourcespath)
          tdfpath <- tempfile(pattern = paste("sim_tdf", i, j, sep = "_"), fileext = ".csv")
          write_csv(sim_food_web$TDF_summary |> 
                      mutate(across(where(is.double), ~ 0)),
                    file = tdfpath)
          consumerspath <- tempfile(pattern = paste("sim_consumers", i, j, sep = "_"), fileext = ".csv")
          write_csv(consumers_rescaled, 
                    file = consumerspath)
          
          mix <- MixSIAR::load_mix_data(consumerspath, 
                                        iso_names = c("d13C", "d15N"), 
                                        factors = NULL,
                                        fac_random = NULL,
                                        fac_nested = NULL,
                                        cont_effects = NULL)
          source <- MixSIAR::load_source_data(sourcespath, 
                                              data_type = "raw", 
                                              mix = mix,
                                              conc_dep = FALSE)
          tdf <- MixSIAR::load_discr_data(tdfpath,
                                          mix = mix)
          
          modelpath <- tempfile(pattern = paste("sim_model", i, j, sep = "_"), fileext = ".txt")
          MixSIAR::write_JAGS_model(modelpath, 
                                    resid_err = ifelse(i == 1, FALSE, TRUE), # process-only when sample size is 1
                                    process_err = TRUE, 
                                    mix = mix, 
                                    source = source)
          
          fit_mixsiar <- run_model_parallel(run = run, 
                                            mix = mix, 
                                            source = source, 
                                            discr = tdf, 
                                            model_filename = modelpath, 
                                            alpha.prior = 1,
                                            seed = j)
          
          # Set output options to avoid saving summary, diagnostics and plots
          mixsiar_options <- list(summary_save = FALSE,
                                  summary_name = "summary_statistics",
                                  sup_post = TRUE,
                                  plot_post_save_pdf = FALSE,
                                  plot_post_name = "posterior_density",
                                  sup_pairs = TRUE,
                                  plot_pairs_save_pdf = FALSE,
                                  plot_pairs_name = "pairs_plot",
                                  sup_xy = TRUE,
                                  plot_xy_save_pdf = FALSE,
                                  plot_xy_name = "xy_plot",
                                  gelman = TRUE,
                                  heidel = FALSE,
                                  geweke = TRUE,
                                  diag_save = FALSE,
                                  diag_name = "diagnostics",
                                  indiv_effect = FALSE,       
                                  plot_post_save_png = FALSE, 
                                  plot_pairs_save_png = FALSE,
                                  plot_xy_save_png = FALSE,
                                  diag_save_ggmcmc = FALSE,
                                  return_obj = TRUE)
          
          # # Get diagnostics
          # diagnostics <- MixSIAR::output_diagnostics(fit_mixsiar, 
          #                                            mix = mix, 
          #                                            source = source, 
          #                                            output_options = mixsiar_options)
          
          # Get statistics
          stats <- MixSIAR::output_stats(fit_mixsiar, 
                                         mix = mix, 
                                         source = source, 
                                         output_options = mixsiar_options) |>
            as.data.frame() |>
            tibble::rownames_to_column("parameter")
          
          out <- stats |> 
            filter(startsWith(parameter, "p.global")) |> 
            tidyr::separate_wider_delim(parameter, delim = ".", names = c("p", "global", "source")) |>
            select(source, alpha_median = `50%`, alpha_lower = `2.5%`, alpha_upper = `97.5%`)
          
          return(out)
        }
      )
    }
  )
  
  # clean tempdir
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".csv"))
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".txt"))
  
  return(sim_sample_size)
  
}

#' Test effect of source sample size on mixing model accuracy
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param sample_size A numeric vector of sample sizes to simulate
#' @param nrep Number of simulations to run
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A list
run_simulations_source_sample_size <- function(sim_food_web,
                                               sample_size = c(1, 3, 5, 10, 20), 
                                               nrep = 100,
                                               run = "test") {
  
  # Estimate TP
  si_data_list <- tRophicPosition::loadIsotopeData(sim_food_web$isotopes |> 
                                                     filter(type != "source"),
                                                   consumer = "family 1_species 1",
                                                   b1 = "primary consumer 1", 
                                                   b2 = "primary consumer 2", 
                                                   baselineColumn = "name", 
                                                   consumersColumn = "name",
                                                   groupsColumn = NULL,
                                                   d13C = "d13C",
                                                   d15N = "d15N",
                                                   deltaC = sim_food_web$TDF$deltaC,
                                                   deltaN = sim_food_web$TDF$deltaN)
  
  model.string <- tRophicPosition::jagsBayesianModel(model = "twoBaselinesFull", lambda = 2)
  
  fit_TP <- tRophicPosition::TPmodel(data = si_data_list, 
                                     model.string = model.string, 
                                     n.adapt = 2000, 
                                     n.chains = 4, 
                                     quiet = FALSE)
  
  post <- tRophicPosition::posteriorTP(fit_TP,
                                       variable.names = c("TP", "alpha"),
                                       n.iter = 2000,
                                       burnin = 2000,
                                       thin = 1,
                                       quiet = FALSE)
  
  post_combined <- bind_rows(purrr::map(post, as.data.frame))
  
  TP <- tRophicPosition::getPosteriorMode(post_combined)$TP
  
  # Rescale consumers
  mean_TDF_C <- mean(sim_food_web$TDF$deltaC)
  mean_TDF_N <- mean(sim_food_web$TDF$deltaN)
  
  consumers_rescaled <- sim_food_web$isotopes |> 
    filter(type == "higher consumer") |> 
    mutate(d13C = d13C - mean_TDF_C * (TP - 1),
           d15N = d15N - mean_TDF_N * (TP - 1)) |> 
    select(d13C, d15N)
  
  # Simulate "nrep" source data sets for a range of sample sizes
  # from source means and SDs
  # Estimate source contribution to consumer and return summary table
  sim_sample_size <- purrr::map(
    sample_size |> 
      rlang::set_names(),
    function(i) {
      purrr::map(
        1:nrep,
        function(j) {
          
          set.seed(j)
          
          # Simulate data
          source_names <- sim_food_web$sources_summary$name
          
          sources <- purrr::map(
            1:length(source_names),
            function(x) {
              data.frame(source = factor(source_names[x]),
                         d13C = rnorm(i, sim_food_web$sources_summary$d13C_mean[x], sim_food_web$sources_summary$d13C_sd[x]),
                         d15N = rnorm(i, sim_food_web$sources_summary$d15N_mean[x], sim_food_web$sources_summary$d15N_sd[x]))
            }
          ) |> 
            bind_rows()
          
          # Mixing model
          sourcespath <- tempfile(pattern = paste("sim_sources", i, j, sep = "_"), fileext = ".csv")
          write_csv(sources, 
                    file = sourcespath)
          tdfpath <- tempfile(pattern = paste("sim_tdf", i, j, sep = "_"), fileext = ".csv")
          write_csv(sim_food_web$TDF_summary |> 
                      mutate(across(where(is.double), ~ 0)),
                    file = tdfpath)
          consumerspath <- tempfile(pattern = paste("sim_consumers", i, j, sep = "_"), fileext = ".csv")
          write_csv(consumers_rescaled, 
                    file = consumerspath)
          
          mix <- MixSIAR::load_mix_data(consumerspath, 
                                        iso_names = c("d13C", "d15N"), 
                                        factors = NULL,
                                        fac_random = NULL,
                                        fac_nested = NULL,
                                        cont_effects = NULL)
          source <- MixSIAR::load_source_data(sourcespath, 
                                              data_type = "raw", 
                                              mix = mix,
                                              conc_dep = FALSE)
          tdf <- MixSIAR::load_discr_data(tdfpath,
                                          mix = mix)
          
          modelpath <- tempfile(pattern = paste("sim_model", i, j, sep = "_"), fileext = ".txt")
          MixSIAR::write_JAGS_model(modelpath, 
                                    resid_err = TRUE,
                                    process_err = TRUE, 
                                    mix = mix, 
                                    source = source)
          
          fit_mixsiar <- run_model_parallel(run = run, 
                                            mix = mix, 
                                            source = source, 
                                            discr = tdf, 
                                            model_filename = modelpath, 
                                            alpha.prior = 1,
                                            seed = j)
          
          # Set output options to avoid saving summary, diagnostics and plots
          mixsiar_options <- list(summary_save = FALSE,
                                  summary_name = "summary_statistics",
                                  sup_post = TRUE,
                                  plot_post_save_pdf = FALSE,
                                  plot_post_name = "posterior_density",
                                  sup_pairs = TRUE,
                                  plot_pairs_save_pdf = FALSE,
                                  plot_pairs_name = "pairs_plot",
                                  sup_xy = TRUE,
                                  plot_xy_save_pdf = FALSE,
                                  plot_xy_name = "xy_plot",
                                  gelman = TRUE,
                                  heidel = FALSE,
                                  geweke = TRUE,
                                  diag_save = FALSE,
                                  diag_name = "diagnostics",
                                  indiv_effect = FALSE,       
                                  plot_post_save_png = FALSE, 
                                  plot_pairs_save_png = FALSE,
                                  plot_xy_save_png = FALSE,
                                  diag_save_ggmcmc = FALSE,
                                  return_obj = TRUE)
          
          # # Get diagnostics
          # diagnostics <- MixSIAR::output_diagnostics(fit_mixsiar, 
          #                                            mix = mix, 
          #                                            source = source, 
          #                                            output_options = mixsiar_options)
          
          # Get statistics
          stats <- MixSIAR::output_stats(fit_mixsiar, 
                                         mix = mix, 
                                         source = source, 
                                         output_options = mixsiar_options) |>
            as.data.frame() |>
            tibble::rownames_to_column("parameter")
          
          out <- stats |> 
            filter(startsWith(parameter, "p.global")) |> 
            tidyr::separate_wider_delim(parameter, delim = ".", names = c("p", "global", "source")) |>
            select(source, alpha_median = `50%`, alpha_lower = `2.5%`, alpha_upper = `97.5%`)
          
          return(out)
        }
      )
    }
  )
  
  # clean tempdir
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".csv"))
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".txt"))
  
  return(sim_sample_size)
  
}
