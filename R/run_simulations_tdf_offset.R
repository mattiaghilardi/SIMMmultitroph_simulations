# Effect of a misspecified mean TDF on estimated source proportions
# (Appendix S4: Section 5, Figure S10).
# In each replicate, the true mean TDF is the mean TDF of `sim_food_web` plus
# an offset drawn from N(0, SD(TDF)). Baselines and consumers are simulated with
# the true TDF; TP is estimated and consumers are rescaled with the assumed
# mean TDF, as in the workflow. Offsets are drawn with seed j + 1e5, so they are
# the same across sample sizes.
# Requires the helpers in run_simulations_tdf_sd.R.

#' Simulate baseline isotope data with a given TDF
#'
#' Same process as the baselines in [simulate_food_web()] (TP = 2, one trophic
#' step above the sources), with the TDF supplied by the user.
#'
#' @param sim_food_web Simulated food web with two baselines
#' @param deltaC,deltaN Numeric vectors of TDFs (their mean and SD are used)
#' @param C Total consumption rate of baselines. Default to 20, as in `_targets.R`
#' @param TP_sd Standard deviation of baseline TP. Default to 0.05
#' @param tdf_variation Passed to [simulate_consumer_data_tdf()]
#'
#' @return A data frame with columns `type`, `name`, `d13C`, `d15N`
simulate_baselines_tdf <- function(sim_food_web, deltaC, deltaN, C = 20, TP_sd = 0.05,
                                   tdf_variation = "none") {

  sim_food_web$alpha_baselines |>
    split(~name, drop = TRUE) |>
    purrr::map(
      ~ simulate_consumer_data_tdf(n = unique(.x$n_obs_baselines),
                                   sources = sim_food_web$sources_summary,
                                   alpha = .x$alpha,
                                   C = C,
                                   TP_mean = 2,
                                   TP_sd = TP_sd,
                                   deltaC = deltaC,
                                   deltaN = deltaN,
                                   tdf_variation = tdf_variation,
                                   seed = NA) # seed set upfront
    ) |>
    bind_rows(.id = "name") |>
    mutate(type = "baseline") |>
    select(type, name, d13C, d15N)
}


#' Test the effect of a misspecified mean TDF on estimated source proportions
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param sample_size A numeric vector of consumer sample sizes. For a sample
#' size of 1, only the process error structure is fitted
#' @param nrep Number of simulations to run
#' @param C Total consumer consumption rate
#' @param C_baselines Total baseline consumption rate (as in `simulate_food_web()`)
#' @param tdf_variation Individual variation in TDF around the true mean:
#' "none", "per_individual" or "per_step" (see [simulate_consumer_data_tdf()])
#' @param offset_sd_C,offset_sd_N SD of the offset in the mean TDF. If NULL
#' (default), the SD of the TDFs in `sim_food_web` (Post 2002)
#' @param tdf_sd_levels TDF SD levels supplied to MixSIAR (see [get_tdf_sd_levels()])
#' @param error_structures Error structures for sample sizes > 1
#' @param TP_prior Prior on consumer TP, passed to [estimate_tp_consumer()]
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A tibble with one row per sample size, replicate, error structure,
#' TDF SD level and source, including the TDF offsets and the true and
#' estimated TP
run_simulations_tdf_offset <- function(sim_food_web,
                                       sample_size = c(1, 5, 10, 20),
                                       nrep = 50,
                                       C = 40,
                                       C_baselines = 20,
                                       tdf_variation = c("none", "per_individual", "per_step"),
                                       offset_sd_C = NULL,
                                       offset_sd_N = NULL,
                                       tdf_sd_levels = c("zero", "one_step", "scaled_TP"),
                                       error_structures = "process_residual",
                                       TP_prior = "dunif(2, 10)",
                                       run = "test") {

  tdf_variation <- rlang::arg_match(tdf_variation)
  all_tdf_sd_levels <- c("zero", "one_step", "scaled_sqrt", "scaled", "scaled_TP")
  tdf_sd_levels <- rlang::arg_match(tdf_sd_levels, all_tdf_sd_levels, multiple = TRUE)
  error_structures <- rlang::arg_match(error_structures,
                                       c("process_residual", "process"),
                                       multiple = TRUE)
  if (any(sample_size < 1)) stop("All `sample_size` values must be >= 1")

  # Assumed (literature) TDFs, used for TP estimation and rescaling
  deltaC <- sim_food_web$TDF$deltaC
  deltaN <- sim_food_web$TDF$deltaN
  if (is.null(offset_sd_C)) offset_sd_C <- sd(deltaC)
  if (is.null(offset_sd_N)) offset_sd_N <- sd(deltaN)

  alpha_true <- sim_food_web$alpha_species |>
    select(source, alpha_true = alpha)
  TP_true <- sim_food_web$TP_species$TP

  # Source data are fixed across replicates
  sourcespath <- tempfile(pattern = "sim_tdf_off_sources", fileext = ".csv")
  write_csv(sim_food_web$isotopes |>
              filter(type == "source") |>
              select(source = name, d13C, d15N),
            file = sourcespath)

  mixsiar_out <- purrr::map(
    sample_size |>
      rlang::set_names(),
    function(i) {
      purrr::map(
        1:nrep,
        function(j) {

          # True mean TDF of the system: shared offset (same across sample sizes)
          set.seed(j + 1e5)
          offset_C <- rnorm(1, 0, offset_sd_C)
          offset_N <- rnorm(1, 0, offset_sd_N)
          deltaC_true <- deltaC + offset_C
          deltaN_true <- deltaN + offset_N

          # Baselines with the true TDF (continue the same random stream)
          baselines <- simulate_baselines_tdf(sim_food_web,
                                              deltaC = deltaC_true,
                                              deltaN = deltaN_true,
                                              C = C_baselines,
                                              tdf_variation = tdf_variation)

          # Consumers with the true TDF
          data <- simulate_consumer_data_tdf(i,
                                             sources = sim_food_web$sources_summary,
                                             alpha = sim_food_web$alpha_species$alpha,
                                             C = C,
                                             TP_mean = TP_true,
                                             TP_sd = 0.05,
                                             deltaC = deltaC_true,
                                             deltaN = deltaN_true,
                                             tdf_variation = tdf_variation,
                                             seed = j)

          # Estimate TP with the assumed TDF
          tp <- estimate_tp_consumer(sim_food_web, data,
                                     TP_prior = TP_prior,
                                     baselines = baselines)

          # Rescale consumers with the assumed mean TDF (as in the workflow)
          consumers_rescaled <- data |>
            mutate(d13C = d13C - mean(deltaC) * (tp$TP - 1),
                   d15N = d15N - mean(deltaN) * (tp$TP - 1))

          # TDF SDs for this replicate (depend on the TP estimate)
          sd_levels <- get_tdf_sd_levels(sim_food_web,
                                         TP = tp$TP,
                                         draws = tp$draws,
                                         levels = tdf_sd_levels)

          # Error structures for this sample size: process only when n = 1
          err_i <- if (i == 1) "process" else error_structures

          fit_mixsiar_tdf_sd(consumers_rescaled,
                             sourcespath = sourcespath,
                             sim_food_web = sim_food_web,
                             sd_levels = sd_levels,
                             error_structures = err_i,
                             id = paste("sim_tdf_off", i, j, sep = "_"),
                             run = run,
                             jags.seed = j) |>
            mutate(offset_C = offset_C,
                   offset_N = offset_N,
                   TP = tp$TP,
                   TP_true = TP_true,
                   .before = 1)
        }
      ) |>
        bind_rows(.id = "replicate")
    }
  ) |>
    bind_rows(.id = "sample_size")

  # clean tempdir
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = "^sim_tdf_off_.*\\.(csv|txt)$"))

  mixsiar_out |>
    left_join(alpha_true, by = "source") |>
    mutate(tdf_variation = tdf_variation,
           sample_size = as.numeric(sample_size),
           replicate = as.integer(replicate),
           TP_error = TP - TP_true,
           # Shift in rescaled consumer values relative to the true mixture:
           # offset * (TP_true - 1) - assumed mean TDF * TP error
           shift_C = offset_C * (TP_true - 1) - mean(deltaC) * TP_error,
           shift_N = offset_N * (TP_true - 1) - mean(deltaN) * TP_error,
           tdf_sd_level = factor(tdf_sd_level, levels = all_tdf_sd_levels),
           error_structure = factor(error_structure,
                                    levels = c("process_residual", "process")),
           .before = 1)
}


#' Expected error in TP from baseline mismatch and TDF offsets
#'
#' Extends the two-baseline derivation of [estimate_tp_error()] (Appendix S7) to
#' a misspecified mean TDF. Expected (population mean) isotope values of the
#' consumer and of the two baselines are computed from the source means with
#' the true TDF (assumed TDF + offset); TP is then solved with the assumed TDF,
#' as in the two-baseline model:
#'   d13C_c = TDF_C * (TP - 2) + a * d13C_b1 + (1 - a) * d13C_b2
#'   d15N_c = TDF_N * (TP - 2) + a * d15N_b1 + (1 - a) * d15N_b2
#' with `a` constrained between 0 and 1. With both offsets equal to 0 this is
#' the baseline-mismatch error of [estimate_tp_error()]. Sampling error is ignored.
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param offset_C,offset_N Numeric vectors (same length) of offsets in the mean TDF
#' @param TP True consumer TP. Default to the species TP in `sim_food_web`
#'
#' @return A numeric vector of expected TP errors (estimated - true TP)
expected_tp_error <- function(sim_food_web, offset_C = 0, offset_N = 0,
                              TP = sim_food_web$TP_species$TP) {
  
  src <- sim_food_web$sources_summary
  TDF_C <- mean(sim_food_web$TDF$deltaC)
  TDF_N <- mean(sim_food_web$TDF$deltaN)
  
  # Source mixtures (no trophic discrimination) of consumer and baselines
  mix <- function(source, alpha, iso) {
    sum(alpha * src[[paste0(iso, "_mean")]][match(as.character(source), 
                                                  as.character(src$name))])
  }
  cons <- sim_food_web$alpha_species
  mix_C_c <- mix(cons$source, cons$alpha, "d13C")
  mix_N_c <- mix(cons$source, cons$alpha, "d15N")
  
  base <- sim_food_web$alpha_baselines |> 
    group_by(name) |> 
    summarise(C = mix(source, alpha, "d13C"),
              N = mix(source, alpha, "d15N"),
              .groups = "drop") |> 
    arrange(name) # primary consumer 1 = b1, primary consumer 2 = b2
  
  # Differences consumer - baseline 2 (true TDF; baselines at TP = 2)
  x <- mix_C_c - base$C[2] + (TDF_C + offset_C) * (TP - 2)
  y <- mix_N_c - base$N[2] + (TDF_N + offset_N) * (TP - 2)
  # Differences baseline 1 - baseline 2 (TDF cancels)
  dC_b <- base$C[1] - base$C[2]
  dN_b <- base$N[1] - base$N[2]
  
  # Solve with the assumed TDF
  r <- TDF_C / TDF_N
  a <- (x - r * y) / (dC_b - r * dN_b)
  a <- pmin(pmax(a, 0), 1)
  TP_hat <- 2 + (y - a * dN_b) / TDF_N
  
  TP_hat - TP
}
