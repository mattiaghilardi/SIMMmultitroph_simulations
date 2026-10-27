# Sensitivity of estimated source proportions to TDF uncertainty
# (Appendix S4: Section 5, Figures S8-S9).
# Consumers are simulated from `sim_food_web` with or without individual
# variation in TDF, TP is estimated and consumers are rescaled as in the
# workflow, and MixSIAR is fitted with different TDF SDs (see
# get_tdf_sd_levels()). Fits within a replicate share the same data, TP
# estimate and JAGS seed.
# Also contains helpers used by run_simulations_tdf_offset.R.

#' Simulate consumer isotope data with individual variation in TDF
#'
#' Same data-generating process as [simulate_consumer_data()], with optional
#' individual variation in the total trophic discrimination. When
#' `tdf_variation = "none"`, it calls [simulate_consumer_data()] directly so that
#' data sets are identical to those of the consumer sample-size simulations.
#'
#' @inheritParams simulate_consumer_data
#' @param tdf_variation One of "none", "per_individual" or "per_step"
#'
#' @return A data frame with `n` rows and two columns, "d13C" and "d15N"
simulate_consumer_data_tdf <- function(n, sources, alpha, C = 10, TP_mean, TP_sd = 0.05,
                                       deltaC, deltaN,
                                       tdf_variation = c("none", "per_individual", "per_step"),
                                       seed = NA) {

  tdf_variation <- rlang::arg_match(tdf_variation)

  if (tdf_variation == "none") {
    return(simulate_consumer_data(n, sources = sources, alpha = alpha, C = C,
                                  TP_mean = TP_mean, TP_sd = TP_sd,
                                  deltaC = deltaC, deltaN = deltaN, seed = seed))
  }

  if (!is.na(seed)) set.seed(seed)

  # Mixture values without trophic discrimination (same process as simulate_consumer_data())
  Zik <- rmultinom(n, C, alpha)

  get_tracers_values <- function(j) {
    j_mean <- sources |> dplyr::pull(paste(j, "mean", sep = "_"))
    j_sd <- sources |> dplyr::pull(paste(j, "sd", sep = "_"))
    Yijk <- purrr::map(
      .x = 1:dim(Zik)[2],
      function(i)
        purrr::map(1:dim(Zik)[1],
                   ~ sum(rnorm(Zik[.x, i], j_mean[.x], j_sd[.x])) / Zik[.x, i]
        ) |>
        purrr::as_vector()
    )
    Xij <- purrr::map_dbl(1:length(Yijk), ~ sum(Yijk[[.x]] * Zik[, .x], na.rm = TRUE) / C)
    res_err <- rnorm(n, 0, 0.1)
    Xij + res_err
  }

  TP_reps <- rnorm(n, TP_mean, TP_sd)
  steps <- TP_reps - 1

  d13C <- get_tracers_values("d13C")
  d15N <- get_tracers_values("d15N")

  # Total trophic discrimination for each individual
  if (tdf_variation == "per_individual") {
    # One TDF per individual, applied at each trophic step
    total_C <- rnorm(n, mean(deltaC), sd(deltaC)) * steps
    total_N <- rnorm(n, mean(deltaN), sd(deltaN)) * steps
  } else {
    # Independent TDF at each trophic step: sum of (TP - 1) independent draws
    total_C <- rnorm(n, mean(deltaC) * steps, sd(deltaC) * sqrt(steps))
    total_N <- rnorm(n, mean(deltaN) * steps, sd(deltaN) * sqrt(steps))
  }

  data.frame(d13C = d13C + total_C, d15N = d15N + total_N)
}


#' Estimate consumer trophic position with tRophicPosition
#'
#' Same model and settings as in [run_simulations_consumer_sample_size()]
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param data Consumer isotope data, output of [simulate_consumer_data_tdf()]
#' @param TP_prior Prior on TP, passed to [tRophicPosition::jagsBayesianModel()].
#' Default to "dunif(2, 10)", the tRophicPosition default used in the other simulations
#' @param baselines Baseline isotope data (columns `type`, `name`, `d13C`, `d15N`).
#' If NULL (default), the baselines of `sim_food_web` are used
#'
#' @return A list with the posterior mode of TP (`TP`) and a data frame of joint
#' posterior draws of TP, muDeltaC and muDeltaN (`draws`)
estimate_tp_consumer <- function(sim_food_web, data, TP_prior = "dunif(2, 10)",
                                 baselines = NULL) {

  if (is.null(baselines)) {
    baselines <- sim_food_web$isotopes |>
      filter(type == "baseline")
  }

  si_data_list <- tRophicPosition::loadIsotopeData(baselines |>
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

  # do.call() inlines the value of TP_prior in the call: tRophicPosition reads
  # its arguments with match.call(), so passing the variable name directly fails
  model.string <- do.call(tRophicPosition::jagsBayesianModel,
                          list(model = "twoBaselinesFull",
                               lambda = 2,
                               TP = TP_prior))

  fit_TP <- tRophicPosition::TPmodel(data = si_data_list,
                                     model.string = model.string,
                                     n.adapt = 2000,
                                     n.chains = 4,
                                     quiet = TRUE)

  post <- tRophicPosition::posteriorTP(fit_TP,
                                       variable.names = c("TP", "alpha",
                                                          "muDeltaC", "muDeltaN"),
                                       n.iter = 2000,
                                       burnin = 2000,
                                       thin = 1,
                                       quiet = TRUE)

  post_combined <- bind_rows(purrr::map(post, as.data.frame))

  list(TP = tRophicPosition::getPosteriorMode(post_combined)$TP,
       draws = post_combined |>
         select(TP, muDeltaC, muDeltaN))
}


#' Compute the TDF SDs supplied to MixSIAR for each level
#'
#' @param sim_food_web Simulated food web
#' @param TP Point estimate (posterior mode) of consumer TP
#' @param draws Joint posterior draws of TP, muDeltaC and muDeltaN, from
#' [estimate_tp_consumer()]
#' @param levels Character vector of TDF SD levels, any of:
#'   - "zero": 0, as in the workflow
#'   - "one_step": SD of the TDF
#'   - "scaled_sqrt": sqrt(TP - 1) * SD
#'   - "scaled": (TP - 1) * SD
#'   - "scaled_TP": sqrt(Var[(TP - 1) * muDelta] + E[(TP - 1)^2] * SD^2), from
#'     joint posterior draws of TP and the mean TDF (Appendix S4: Equation S3)
#'
#' @return A tibble with columns `tdf_sd_level`, `sd_C`, `sd_N`
get_tdf_sd_levels <- function(sim_food_web, TP, draws, levels) {

  sd_C <- sd(sim_food_web$TDF$deltaC)
  sd_N <- sd(sim_food_web$TDF$deltaN)

  # Total SD from joint posterior draws:
  # Var[(TP - 1) * muDelta] + E[(TP - 1)^2] * Var(TDF)
  total_sd <- function(mu_delta, sd_delta) {
    steps <- draws$TP - 1
    sqrt(var(steps * mu_delta) + mean(steps^2) * sd_delta^2)
  }

  out <- tibble::tibble(tdf_sd_level = levels) |>
    mutate(sd_C = NA_real_, sd_N = NA_real_)

  for (l in levels) {
    sds <- switch(
      l,
      zero = c(0, 0),
      one_step = c(sd_C, sd_N),
      scaled_sqrt = sqrt(TP - 1) * c(sd_C, sd_N),
      scaled = (TP - 1) * c(sd_C, sd_N),
      scaled_TP = c(total_sd(draws$muDeltaC, sd_C),
                    total_sd(draws$muDeltaN, sd_N)),
      stop("Unknown TDF SD level: ", l)
    )
    out$sd_C[out$tdf_sd_level == l] <- sds[1]
    out$sd_N[out$tdf_sd_level == l] <- sds[2]
  }

  out
}


#' Extract summaries of source proportions and multiplicative residual error terms
#'
#' @param fit Output of [run_model_parallel()]
#' @param source Output of [MixSIAR::load_source_data()]
#'
#' @return A tibble with one row per source
extract_mixsiar_summary <- function(fit, source) {

  sims <- fit$BUGSoutput$sims.list

  p <- sims$p.global
  colnames(p) <- source$source_names

  out <- tibble::tibble(source = colnames(p),
                        alpha_median = apply(p, 2, median),
                        alpha_lower = apply(p, 2, quantile, probs = 0.025),
                        alpha_upper = apply(p, 2, quantile, probs = 0.975))

  # Multiplicative residual error terms (Epsilon), only for process x residual
  # Columns follow `iso_names` in load_mix_data(): 1 = d13C, 2 = d15N
  if (!is.null(sims$resid.prop)) {
    eps <- sims$resid.prop
    out <- out |>
      mutate(eps_C_median = median(eps[, 1]),
             eps_C_lower = quantile(eps[, 1], 0.025),
             eps_C_upper = quantile(eps[, 1], 0.975),
             eps_N_median = median(eps[, 2]),
             eps_N_lower = quantile(eps[, 2], 0.025),
             eps_N_upper = quantile(eps[, 2], 0.975))
  } else {
    out <- out |>
      mutate(eps_C_median = NA_real_, eps_C_lower = NA_real_, eps_C_upper = NA_real_,
             eps_N_median = NA_real_, eps_N_lower = NA_real_, eps_N_upper = NA_real_)
  }

  out
}


#' Fit MixSIAR to rescaled consumers for each error structure and TDF SD level
#'
#' Shared by [run_simulations_tdf_sd()] and [run_simulations_tdf_offset()].
#' TDF means are set to 0 (consumers already rescaled) and TDF SDs to the
#' values in `sd_levels`. All fits use the same JAGS seed.
#'
#' @param consumers_rescaled Data frame with rescaled consumer `d13C` and `d15N`
#' @param sourcespath Path to the csv file with raw source data
#' @param sim_food_web Simulated food web (for the TDF file template)
#' @param sd_levels Output of [get_tdf_sd_levels()]
#' @param error_structures Error structures to fit: "process_residual" and/or "process"
#' @param id A string used to name temporary files
#' @param run The `run` argument passed to [run_model_parallel()]
#' @param jags.seed Seed passed to [run_model_parallel()]
#'
#' @return A tibble with one row per error structure, TDF SD level and source
fit_mixsiar_tdf_sd <- function(consumers_rescaled, sourcespath, sim_food_web,
                               sd_levels, error_structures, id, run, jags.seed) {

  consumerspath <- tempfile(pattern = paste(id, "consumers", sep = "_"),
                            fileext = ".csv")
  write_csv(consumers_rescaled |> select(d13C, d15N), file = consumerspath)

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

  purrr::map(
    error_structures |>
      rlang::set_names(),
    function(err) {

      modelpath <- tempfile(pattern = paste(id, "model", err, sep = "_"),
                            fileext = ".txt")
      MixSIAR::write_JAGS_model(modelpath,
                                resid_err = err == "process_residual",
                                process_err = TRUE,
                                mix = mix,
                                source = source)

      purrr::map(
        seq_len(nrow(sd_levels)),
        function(l) {

          # TDF means = 0 (consumers already rescaled), SDs = level value
          tdfpath <- tempfile(pattern = paste(id, "tdf", err, l, sep = "_"),
                              fileext = ".csv")
          write_csv(sim_food_web$TDF_summary |>
                      mutate(Meand13C = 0,
                             SDd13C = sd_levels$sd_C[l],
                             Meand15N = 0,
                             SDd15N = sd_levels$sd_N[l]),
                    file = tdfpath)
          tdf <- MixSIAR::load_discr_data(tdfpath, mix = mix)

          fit_mixsiar <- run_model_parallel(run = run,
                                            mix = mix,
                                            source = source,
                                            discr = tdf,
                                            model_filename = modelpath,
                                            alpha.prior = 1,
                                            jags.seed = jags.seed)

          extract_mixsiar_summary(fit_mixsiar, source) |>
            mutate(tdf_sd_level = sd_levels$tdf_sd_level[l],
                   sd_C = sd_levels$sd_C[l],
                   sd_N = sd_levels$sd_N[l],
                   .before = 1)
        }
      ) |>
        bind_rows()
    }
  ) |>
    bind_rows(.id = "error_structure")
}


#' Test the effect of the TDF SD supplied to MixSIAR on estimated source proportions
#'
#' @param sim_food_web Simulated food web with two baselines and one consumer species
#' @param sample_size A numeric vector of consumer sample sizes. For a sample
#' size of 1, only the process error structure is fitted
#' @param nrep Number of simulations to run
#' @param C Total consumer consumption rate
#' @param tdf_variation Individual variation in TDF in the simulated consumers:
#' "none", "per_individual" or "per_step"
#' @param tdf_sd_levels TDF SD levels supplied to MixSIAR, any of "zero",
#' "one_step", "scaled_sqrt", "scaled", "scaled_TP" (see [get_tdf_sd_levels()]). Default
#' to all five; choose a subset in the call
#' @param error_structures MixSIAR error structures: "process_residual" and/or "process"
#' @param TP_prior Prior on consumer TP, passed to [estimate_tp_consumer()].
#' At n = 1 the TP posterior, and hence the "scaled_TP" level, depends strongly
#' on this prior
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A tibble with one row per sample size, replicate, error structure,
#' TDF SD level and source
run_simulations_tdf_sd <- function(sim_food_web,
                                   sample_size = 10,
                                   nrep = 50,
                                   C = 40,
                                   tdf_variation = c("none", "per_individual", "per_step"),
                                   tdf_sd_levels = c("zero", "one_step", "scaled_sqrt",
                                                     "scaled", "scaled_TP"),
                                   error_structures = c("process_residual", "process"),
                                   TP_prior = "dunif(2, 10)",
                                   run = "test") {

  tdf_variation <- rlang::arg_match(tdf_variation)
  all_tdf_sd_levels <- c("zero", "one_step", "scaled_sqrt", "scaled", "scaled_TP")
  tdf_sd_levels <- rlang::arg_match(tdf_sd_levels, all_tdf_sd_levels, multiple = TRUE)
  error_structures <- rlang::arg_match(error_structures,
                                       c("process_residual", "process"),
                                       multiple = TRUE)
  if (any(sample_size < 1)) stop("All `sample_size` values must be >= 1")
  if (any(sample_size == 1) && !"process" %in% error_structures) {
    stop("`sample_size = 1` requires the \"process\" error structure ",
         "(the residual error term cannot be estimated from a single observation)")
  }
  if (any(sample_size == 1) && "process_residual" %in% error_structures) {
    message("Only process error is fitted for `sample_size = 1`")
  }

  alpha_true <- sim_food_web$alpha_species |>
    select(source, alpha_true = alpha)

  # Source data are fixed across replicates
  sourcespath <- tempfile(pattern = "sim_tdf_sd_sources", fileext = ".csv")
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

          # Simulate data
          data <- simulate_consumer_data_tdf(i,
                                             sources = sim_food_web$sources_summary,
                                             alpha = sim_food_web$alpha_species$alpha,
                                             C = C,
                                             TP_mean = sim_food_web$TP_species$TP,
                                             TP_sd = 0.05,
                                             deltaC = sim_food_web$TDF$deltaC,
                                             deltaN = sim_food_web$TDF$deltaN,
                                             tdf_variation = tdf_variation,
                                             seed = j)

          # Estimate TP (once per replicate)
          tp <- estimate_tp_consumer(sim_food_web, data, TP_prior = TP_prior)

          # Rescale consumers with the mean TDF (as in the workflow)
          mean_TDF_C <- mean(sim_food_web$TDF$deltaC)
          mean_TDF_N <- mean(sim_food_web$TDF$deltaN)

          consumers_rescaled <- data |>
            mutate(d13C = d13C - mean_TDF_C * (tp$TP - 1),
                   d15N = d15N - mean_TDF_N * (tp$TP - 1))

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
                             id = paste("sim_tdf_sd", i, j, sep = "_"),
                             run = run,
                             jags.seed = j) |>
            mutate(TP = tp$TP,
                   TP_true = sim_food_web$TP_species$TP,
                   .before = 1)
        }
      ) |>
        bind_rows(.id = "replicate")
    }
  ) |>
    bind_rows(.id = "sample_size")

  # clean tempdir
  file.remove(list.files(tempdir(), full.names = TRUE, pattern = "^sim_tdf_sd_.*\\.(csv|txt)$"))

  mixsiar_out |>
    left_join(alpha_true, by = "source") |>
    mutate(tdf_variation = tdf_variation,
           sample_size = as.numeric(sample_size),
           replicate = as.integer(replicate),
           TP_error = TP - TP_true,
           # fixed levels (canonical order) so outputs with different subsets combine
           tdf_sd_level = factor(tdf_sd_level, levels = all_tdf_sd_levels),
           # fixed levels so outputs from separate runs combine consistently
           error_structure = factor(error_structure,
                                    levels = c("process_residual", "process")),
           .before = 1)
}
