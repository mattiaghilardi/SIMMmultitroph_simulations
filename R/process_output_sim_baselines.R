#' Compute the absolute error in *p~k~* and the proportion of true *p~k~* within the 95% CI
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param nrep Number of simulations run
#'
#' @return A data frame
get_performance_pk <- function(sim_baselines, nrep = 100) {
  
  purrr::map(
  1:2,
  function(i)
    purrr::map(1:nrep,
               function(j)
                 sim_baselines[[i]][[j]][["estimated_pk"]] |>
                 mutate(n_baselines = i, replicate = j))) |>
  bind_rows() |>
  mutate(n_baselines = as.factor(n_baselines),
         within_95CI = ifelse(alpha_observed >= alpha_lower & alpha_observed <= alpha_upper, 1, 0),
         abs_error = abs(alpha_median - alpha_observed))
}

#' Compute the absolute error in TP and the proportion of true TP within the 95% CI
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param nrep Number of simulations run
#'
#' @return A data frame
get_performance_TP <- function(sim_baselines, nrep = 100) {
  
  purrr::map(
  1:2,
  function(i)
    purrr::map(1:nrep,
               function(j)
                 sim_baselines[[i]][[j]][["estimated_TP"]] |>
                 mutate(n_baselines = i, replicate = j))) |>
  bind_rows() |>
  mutate(n_baselines = as.factor(n_baselines),
         within_95CI = ifelse(TP_observed >= lower & TP_observed <= upper, 1, 0),
         abs_error = abs(mode - TP_observed))
}

#' Join absolute TP error to the performance of source proportions
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#'
#' @return `performance_pk` with an `abs_error_TP` column
join_tp_error <- function(performance_pk, performance_TP) {
  performance_pk |>
    left_join(performance_TP |>
                select(n_baselines, replicate, abs_error_TP = abs_error),
              by = c("n_baselines", "replicate"))
}

#' Maximum absolute error in *p~k~* per food web, with absolute TP error
#'
#' @inheritParams join_tp_error
#'
#' @return A tibble with one row per number of baselines and replicate
get_max_abs_error <- function(performance_pk, performance_TP) {
  join_tp_error(performance_pk, performance_TP) |>
    group_by(n_baselines, replicate, abs_error_TP) |>
    summarise(max_abs_error = max(abs_error)) |>
    ungroup()
}

#' Model the maximum absolute error in *p~k~*
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param seed Seed for reproducibility
#'
#' @return A `brmsfit` object
fit_baseline_effect <- function(performance_pk, performance_TP, seed = NA) {
  
  data <- get_max_abs_error(performance_pk, performance_TP) |> 
    mutate(replicate = as.factor(replicate))
  
  fit <- brm(bf(max_abs_error ~ 1 + n_baselines * log(abs_error_TP) + (1 | replicate),
                phi ~ 1 + n_baselines + (1 | replicate)),
             family = Beta(),
             data = data,
             prior = prior(normal(0, 1), class = b),
             cores = 4,
             backend = "cmdstan",
             control = list(adapt_delta = 0.999),
             seed = seed)
  
  plot(fit, ask = FALSE)
  
  return(fit)
}

#' Model coverage of 95% CI as a function of absolute TP error
#'
#' Logistic regression of whether each interval contains the true proportion,
#' with the same structure as [fit_baseline_effect()]: interaction between the
#' number of baselines and log absolute TP error, and a random intercept per
#' replicate (the same food web is simulated with one and two baselines)
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param seed Seed for reproducibility
#'
#' @return A `brmsfit` object
fit_coverage_tp_error <- function(performance_pk, performance_TP, seed = NA) {
  
  data <- join_tp_error(performance_pk, performance_TP) |>
    mutate(covered = as.integer(within_95CI),
           n_baselines = factor(n_baselines),
           replicate = as.factor(replicate))
  
  fit <- brm(covered ~ 1 + n_baselines * log(abs_error_TP) + (1 | replicate),
             family = bernoulli(),
             data = data,
             prior = prior(normal(0, 1), class = b),
             cores = 4,
             backend = "cmdstan",
             seed = seed)
  
  plot(fit, ask = FALSE)
  
  return(fit)
}

#' Extract contrasts between one and two baselines in model parameters
#'
#' Posterior contrasts (one baseline - two baselines) for the slope of log
#' absolute TP error, in both the model of maximum absolute error in *p~k~* and
#' the model of coverage, and for the precision (phi) of maximum absolute error
#' (which does not depend on TP error). Contrasts in mean and coverage depend on
#' TP error because of the interaction and are computed with
#' [get_baseline_contrast_curves()]
#'
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param model_coverage Output of [fit_coverage_tp_error()]
#' @param seed Seed for reproducibility
#'
#' @return A list of two [brms::brmshypothesis] objects: `max_error` (slope and
#' precision) and `coverage` (slope)
get_baseline_contrasts <- function(model_baseline, model_coverage, seed = NULL) {
  
  slope <- "logabs_error_TP > logabs_error_TP + n_baselines2:logabs_error_TP"
  
  list(
    max_error = brms::hypothesis(
      model_baseline,
      c("slope" = slope,
        "phi" = "phi_Intercept > phi_Intercept + phi_n_baselines2"),
      seed = seed),
    coverage = brms::hypothesis(
      model_coverage,
      c("slope" = slope),
      seed = seed)
  )
}

#' Slopes of log absolute TP error for one and two baselines
#'
#' @param model A `brmsfit` object: output of [fit_baseline_effect()] or
#' [fit_coverage_tp_error()]
#' @param .width Width of the credible interval
#'
#' @return A tibble with the median and CI of the slope (link scale) for each
#' number of baselines
get_baseline_slopes <- function(model, .width = 0.95) {
  
  brms::as_draws_df(model) |>
    tibble::as_tibble() |>
    dplyr::transmute(`1` = b_logabs_error_TP,
                     `2` = b_logabs_error_TP + `b_n_baselines2:logabs_error_TP`) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "n_baselines",
                        values_to = "slope") |>
    dplyr::group_by(n_baselines) |>
    ggdist::median_qi(slope, .width = .width)
}

#' Contrasts between one and two baselines across absolute TP error
#'
#' Expected maximum absolute error in *p~k~* (mean of the Beta distribution) and
#' expected coverage of the 95% CI of *p~k~* for a typical food web (random
#' effect set to 0), contrasted between one and two baselines along a grid of
#' absolute TP error
#'
#' @param performance_TP Output of [get_performance_TP()]
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param model_coverage Output of [fit_coverage_tp_error()]
#' @param seed Seed for reproducibility
#'
#' @return A tibble of posterior draws of the contrast (`.epred`, one baseline -
#' two baselines) for each response (`var`: "max_error" or "coverage") and value
#' of absolute TP error
get_baseline_contrast_curves <- function(performance_TP,
                                         model_baseline,
                                         model_coverage,
                                         seed = NULL) {
  
  nd <- expand.grid(n_baselines = factor(1:2),
                    abs_error_TP = exp(seq(log(min(performance_TP$abs_error)),
                                           log(max(performance_TP$abs_error)),
                                           length.out = 51)))
  
  list(max_error = model_baseline, coverage = model_coverage) |>
    purrr::map(
      ~ tidybayes::add_epred_draws(nd, .x, re_formula = NA, seed = seed) |>
        tidybayes::compare_levels(.epred,
                                  by = n_baselines,
                                  comparison = list(c("1", "2"))) |>  # one - two
        dplyr::ungroup()
    ) |>
    dplyr::bind_rows(.id = "var")
}

#' Plot contrasts between one and two baselines (Appendix S4: Figure S5)
#'
#' @param baseline_contrast_curves Output of [get_baseline_contrast_curves()]
#' @param baseline_contrasts Output of [get_baseline_contrasts()]
#' @param .width The `.width` argument passed to [ggdist::stat_lineribbon()]
#'
#' @return A `patchwork` object with A) the contrast in mean maximum absolute
#' error in *p~k~* across absolute TP error, B) the contrast in precision of
#' maximum absolute error in *p~k~*, and C) the contrast in coverage of the 95%
#' CI of *p~k~* across absolute TP error
plot_baseline_contrasts <- function(baseline_contrast_curves,
                                    baseline_contrasts,
                                    .width = ppoints(20)) {
  
  # Contrast across absolute TP error
  plot_curve <- function(data, ylab) {
    ggplot(data, aes(x = abs_error_TP, y = .epred)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      ggdist::stat_lineribbon(aes(fill_ramp = after_stat(.width)),
                              .width = .width,
                              alpha = 0.3,
                              fill = "#2171b5") +
      ggdist::scale_fill_ramp_continuous(range = c(1, 0),
                                         guide = ggdist::guide_rampbar(title = "CI",
                                                                       to = "#2171b5")) +
      labs(x = "Absolute error in TP<i><sub>c</sub></i>",
           y = ylab) +
      theme_bw() +
      theme(axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown())
  }
  
  pA <- baseline_contrast_curves |>
    dplyr::filter(var == "max_error") |>
    plot_curve("Contrast in mean max<br>absolute error in <i>p<sub>k</sub></i>") +
    theme(axis.title.x = element_blank())
  
  pC <- baseline_contrast_curves |>
    dplyr::filter(var == "coverage") |>
    plot_curve("Contrast in coverage<br>of 95% CI of <i>p<sub>k</sub></i>")
  
  # Contrast in precision (independent of TP error)
  hyp <- baseline_contrasts$max_error
  phi <- hyp$samples[[which(hyp$hypothesis$Hypothesis == "phi")]]
  
  pB <- tibble::tibble(value = phi) |>
    ggplot(aes(y = value)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggdist::stat_slabinterval(.width = c(0.5, 0.95),
                              slab_alpha = 0.5) +
    labs(x = NULL,
         y = "Log ratio in precision of max absolute error in <i>p<sub>k</sub></i>") +
    theme_bw() +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  design <- "
  AB
  CB
  "
  
  pA + pB + pC +
    patchwork::plot_layout(design = design,
                           widths = c(3, 1),
                           guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")
}

#' Summarise the contrast in coverage between one and two baselines
#'
#' @param baseline_contrast_curves Output of [get_baseline_contrast_curves()]
#' @param .width Width of the credible interval
#'
#' @return A list with `peak` (largest median contrast, its CI and absolute TP
#' error) and `onset` (smallest absolute TP error above which the CI of the
#' contrast excludes zero over the rest of the range)
summarise_coverage_contrast <- function(baseline_contrast_curves, .width = 0.95) {
  
  summ <- baseline_contrast_curves |>
    dplyr::filter(var == "coverage") |>
    dplyr::group_by(abs_error_TP) |>
    ggdist::median_qi(.epred, .width = .width)
  
  # credible from here onwards: lower bound > 0 at this and all larger errors
  credible_onwards <- rev(cumprod(rev(summ$.lower > 0))) == 1
  
  list(
    peak  = dplyr::slice_max(summ, .epred, n = 1),
    onset = summ$abs_error_TP[which(credible_onwards)[1]]
  )
}

#' Estimate TP error
#'
#' @param sim_results Output of [run_simulations_baselines()] for a single simulation
#'
#' @return A value
get_tp_error <- function(sim_results) {
  
  # sources summary
  sources <- sim_results$data$sources_summary
  # baselines summary
  baselines <- sim_results$data$isotopes |> 
    filter(type == "baseline") |> 
    group_by(name) |> 
    summarise(d13C_mean = mean(d13C),
              d15N_mean = mean(d15N))
  # mean TDFs
  TDF_C <- mean(sim_results$data$TDF$deltaC)
  TDF_N <- mean(sim_results$data$TDF$deltaN)
  # source relative contributions to consumer
  consumer_p <- sim_results$data$alpha_species$alpha
  
  estimate_tp_error(sources, baselines, TP_b = 2, TDF_N, TDF_C, consumer_p)
  
}

#' Coverage of 95% CI of source proportions in bins of absolute TP error
#'
#' Replicates are split into `n_bins` bins with (approximately) equal numbers of
#' replicates, based on quantiles of absolute TP error, separately for one and
#' two baselines.
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param n_bins Number of bins
#'
#' @return A tibble with one row per number of baselines and bin, with the mean,
#' minimum and maximum absolute TP error in the bin, coverage, its 95% binomial
#' CI, number of intervals and number of replicates
get_coverage_tp_error <- function(performance_pk, performance_TP, n_bins = 5) {
  
  # Bins of equal numbers of replicates (one TP error per replicate)
  bins <- performance_TP |>
    select(n_baselines, replicate, abs_error_TP = abs_error) |>
    group_by(n_baselines) |>
    mutate(bin = ntile(abs_error_TP, n_bins)) |>
    ungroup() |>
    select(n_baselines, replicate, bin)
  
  join_tp_error(performance_pk, performance_TP) |>
    left_join(bins, by = c("n_baselines", "replicate")) |>
    group_by(n_baselines, bin) |>
    summarise(abs_error_TP_mean = mean(abs_error_TP),
              abs_error_TP_min = min(abs_error_TP),
              abs_error_TP_max = max(abs_error_TP),
              k = sum(within_95CI),
              n = n(),
              nrep = n_distinct(replicate),
              .groups = "drop") |>
    mutate(coverage = k / n,
           # binomial CI treats intervals as independent; the 3 per
           # replicate are not, so these are slightly too narrow
           lower = purrr::map2_dbl(k, n, ~ binom.test(.x, .y)$conf.int[1]),
           upper = purrr::map2_dbl(k, n, ~ binom.test(.x, .y)$conf.int[2]))
}

#' Plot the effect of absolute error in TP on accuracy and coverage
#'
#' @param performance_pk Output of [get_performance_pk()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param model_baseline Output of [fit_baseline_effect()]
#' @param coverage_tp_error Output of [get_coverage_tp_error()]
#' @param model_coverage Output of [fit_coverage_tp_error()]
#' @param .width The `.width` argument passed to [ggdist::stat_lineribbon()]
#' @param alpha_points Point transparency
#' @param alpha_rug Rug transparency
#' @param col_bins Colour of binned coverage and its intervals
#' @param seed Seed for reproducibility
#'
#' @return A `patchwork` object with (A) the maximum absolute error in source
#' proportions and (B) the coverage of their 95% CI against absolute TP error.
#' Panel B shows the fitted coverage for a typical food web (curve and ribbon;
#' random effect set to 0), binned coverage (points: mean TP error of the bin;
#' vertical bars: 95% binomial CI; horizontal bars: range of TP error in the
#' bin) and the raw data as a rug (top: intervals including the true
#' proportion; bottom: intervals excluding it)
plot_TP_error_effect <- function(performance_pk, 
                                 performance_TP, 
                                 model_baseline,
                                 coverage_tp_error,
                                 model_coverage,
                                 .width = ppoints(20), 
                                 alpha_points = 0.2, 
                                 alpha_rug = 0.2,
                                 col_bins = "grey50",
                                 seed = NULL) {
  
  # A: maximum absolute error in pk
  # Fitted values
  preds <- expand.grid(n_baselines = 1:2, 
                       abs_error_TP = exp(seq(log(min(performance_TP$abs_error)),
                                              log(max(performance_TP$abs_error)),
                                              length.out = 51))) |>
    tidybayes::add_epred_draws(model_baseline,
                               re_formula = NA,
                               seed = seed) |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines")))
  
  # Plot
  p1 <- performance_pk |> 
    group_by(n_baselines, replicate) |> 
    summarise(max_abs_error = max(abs_error)) |>
    ungroup() |> 
    left_join(performance_TP |> 
                select(n_baselines, replicate, abs_error_TP = abs_error)) |> 
    mutate(replicate = as.factor(replicate),
           n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines"))) |> 
    ggplot(aes(x = abs_error_TP, y = max_abs_error)) + 
    ggdist::stat_lineribbon(data = preds,
                            aes(y = .epred, 
                                fill_ramp = after_stat(.width)),
                            .width = .width,
                            alpha = 0.3,
                            fill = "#2171b5") +
    ggdist::scale_fill_ramp_continuous(range = c(1, 0), 
                                       guide = ggdist::guide_rampbar(title = "CI",
                                                                     to = "#2171b5")) +
    geom_point(alpha = alpha_points) + 
    facet_grid(~ n_baselines, scales = "free") +
    labs(x = "Absolute error in TP<i><sub>c</sub></i>",
         y = "Max absolute error in <i>p<sub>k</sub></i>") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  # B: coverage of 95% CI of pk
  # Raw data: one row per interval (source x replicate), with absolute TP error
  raw <- join_tp_error(performance_pk, performance_TP) |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines")))
  
  # Binned coverage (bins with equal numbers of replicates)
  bins <- coverage_tp_error |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines")))
  
  # Fitted values
  preds_cov <- expand.grid(n_baselines = factor(1:2),
                           abs_error_TP = exp(seq(log(min(performance_TP$abs_error)),
                                                  log(max(performance_TP$abs_error)),
                                                  length.out = 51))) |>
    tidybayes::add_epred_draws(model_coverage,
                               re_formula = NA,
                               seed = seed) |>
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines")))
  
  # Plot
  p2 <- bins |> 
    ggplot(aes(x = abs_error_TP_mean, y = coverage)) +
    ggdist::stat_lineribbon(data = preds_cov,
                            aes(x = abs_error_TP, y = .epred,
                                fill_ramp = after_stat(.width)),
                            .width = .width,
                            alpha = 0.3,
                            fill = "#2171b5") +
    ggdist::scale_fill_ramp_continuous(range = c(1, 0),
                                       guide = ggdist::guide_rampbar(title = "CI",
                                                                     to = "#2171b5")) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    # raw data: covered intervals on top, not covered at the bottom
    geom_rug(data = raw |> filter(within_95CI == 1),
             aes(x = abs_error_TP), 
             inherit.aes = FALSE,
             sides = "t", 
             alpha = alpha_rug) +
    geom_rug(data = raw |> filter(within_95CI == 0),
             aes(x = abs_error_TP), 
             inherit.aes = FALSE,
             sides = "b", 
             alpha = alpha_rug) +
    # binned coverage: mean TP error of the bin, 95% binomial CI, range of TP error
    geom_linerange(aes(xmin = abs_error_TP_min, xmax = abs_error_TP_max),
                   color = col_bins) +
    geom_linerange(aes(ymin = lower, ymax = upper),
                   color = col_bins) +
    geom_point(shape = 21, fill = col_bins, color = "white",
               size = 2.5) +
    facet_grid(~ n_baselines, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Absolute error in TP<i><sub>c</sub></i>",
         y = "Coverage of 95% CI of <i>p<sub>k</sub></i>") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  # Combine plots
  (p1 / p2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A")
}


#' Plot the observed TP error against the analytical estimate
#'
#' Error in estimated TP against the TP-corrected difference in d15N between
#' consumer and baseline(s), computed with [get_tp_error()].
#'
#' @param sim_baselines Output of [run_simulations_baselines()]
#' @param performance_TP Output of [get_performance_TP()]
#' @param alpha_points Point transparency
#'
#' @return A ggplot object
plot_TP_error_analytic <- function(sim_baselines,
                                   performance_TP,
                                   alpha_points = 0.2) {
  
  # Get TP error
  tp_error <- purrr::map(
    1:2,
    function(i) {
      purrr::map(
        1:100,
        function(j) {
          data.frame(n_baselines = as.factor(i), 
                     replicate = j,
                     estimated_tp_error = get_tp_error(sim_baselines[[i]][[j]]))
        }
      )
    }
  ) |>
    bind_rows()
  
  # Plot relationship between error in TP and difference in d15N
  performance_TP |> 
    mutate(tp_error = mode - TP_observed) |> 
    select(n_baselines, replicate, tp_error) |> 
    left_join(tp_error) |> 
    mutate(n_baselines = factor(n_baselines,
                                levels = c(1, 2),
                                labels = c("One baseline", "Two baselines"))) |> 
    ggplot(aes(x = estimated_tp_error, y = tp_error)) + 
    geom_abline() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = alpha_points) + 
    facet_grid(~ n_baselines, scales = "free") +
    labs(x = "TP-corrected &Delta;&delta;<sup>15</sup>N between consumer and baseline(s)",
         y = "Error in TP<i><sub>c</sub></i>") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
}