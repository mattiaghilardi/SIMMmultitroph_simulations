#' Test accuracy in estimation of source contribution using one or two baselines
#'
#' @param nrep Number of simulations to run
#' @param C_baselines,C_consumers Total consumption rate of baselines and consumers
#' (i.e. biomass taken from all sources in the tissue turnover period). Default to 20
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A list
run_simulations_baselines <- function(nrep = 100, 
                                      C_baselines = 20,
                                      C_consumers = 20,
                                      run = "test") {
  
  sim_baselines <- purrr::map(
    1:2 |> 
      rlang::set_names(),
    function(i) {
      purrr::map(
        1:nrep,
        function(j) {
          # Simulate data
          set.seed(j) # Ensure source SDs are identical for 1 and 2 baselines
          data <- simulate_food_web(n_sources = 3,
                                    n_baselines = i,
                                    n_families = 1,
                                    n_species = 1,
                                    source_C_sd = runif(3, 0.5, 1.5),
                                    source_N_sd = runif(3, 0.5, 1.5),
                                    C_baselines = C_baselines,
                                    C_consumers = C_consumers,
                                    seed = j)
          
          # Estimate TP
          baseline2 <- if (i == 2) "primary consumer 2" else NULL
          model_name <- if (i == 1) "oneBaseline" else "twoBaselinesFull"
          si_data_list <- tRophicPosition::extractIsotopeData(data$isotopes |> 
                                                                filter(type != "source"), 
                                                              b1 = "primary consumer 1", 
                                                              b2 = baseline2,
                                                              baselineColumn = "name", 
                                                              consumersColumn = "name",
                                                              groupsColumn = NULL,
                                                              d13C = "d13C", 
                                                              d15N = "d15N",
                                                              deltaC = data$TDF$deltaC,
                                                              deltaN = data$TDF$deltaN)
          
          estimated_TP <- tRophicPosition::multiSpeciesTP(siDataList = si_data_list,
                                                          lambda = 2,
                                                          model = model_name,
                                                          n.adapt = 2000, n.iter = 2000,
                                                          burnin = 2000, thin = 1, 
                                                          n.chains = 4, quiet = FALSE)
          
          # Rescale consumers
          mean_TDF_C <- mean(data$TDF$deltaC)
          mean_TDF_N <- mean(data$TDF$deltaN)
          
          consumers_rescaled <- data$isotopes |> 
            filter(type == "higher consumer") |> 
            left_join(estimated_TP$df |>
                        select(consumer, TP = mode), 
                      by = c("name" = "consumer")) |> 
            mutate(d13C = d13C - mean_TDF_C * (TP - 1),
                   d15N = d15N - mean_TDF_N * (TP - 1)) |> 
            select(family, species, d13C, d15N, TP)
          
          # Mixing model
          sourcespath <- tempfile(pattern = "sim_sources", fileext = ".csv")
          write_csv(data$isotopes |> 
                      filter(type == "source") |> 
                      select(source = name, d13C, d15N), 
                    file = sourcespath)
          tdfpath <- tempfile(pattern = "sim_tdf", fileext = ".csv")
          write_csv(data$TDF_summary |> 
                      mutate(across(where(is.double), ~ 0)),
                    file = tdfpath)
          consumerspath <- tempfile(pattern = "sim_consumers", fileext = ".csv")
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
          
          modelpath <- tempfile(pattern = "sim_model", fileext = ".txt")
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
          
          # Get diagnostics
          diagnostics <- MixSIAR::output_diagnostics(fit_mixsiar, 
                                                     mix = mix, 
                                                     source = source, 
                                                     output_options = mixsiar_options)
          
          # Get statistics
          stats <- MixSIAR::output_stats(fit_mixsiar, 
                                         mix = mix, 
                                         source = source, 
                                         output_options = mixsiar_options) |>
            as.data.frame() |>
            tibble::rownames_to_column("parameter")
          
          pk <- stats |> 
            filter(startsWith(parameter, "p.global")) |> 
            tidyr::separate_wider_delim(parameter, delim = ".", names = c("p", "global", "source")) |>
            left_join(data$alpha_species) |> 
            select(source, alpha_median = `50%`, alpha_lower = `2.5%`, alpha_upper = `97.5%`, alpha_observed = alpha)
          
          TP <- estimated_TP$df |> mutate(TP_observed = data$TP_species$TP)
          
          return(list(data = data,
                      model_diagnostics = diagnostics,
                      estimated_pk = pk, 
                      estimated_TP = TP))
        }
      )
    }
  )
  
  # clean tempdir
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".csv"))
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = ".txt"))
  
  return(sim_baselines)
}